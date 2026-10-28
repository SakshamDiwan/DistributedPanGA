#!/usr/bin/env python3
"""Alignment-stage balance report for one pga-mpi-merge run directory.

Reads the "[pr1] rank r:" lines from <dir>/merge.err and the per-task records
in <dir>/pairs.rank*.tsv (one row per align_contigs call: type, anchor,
partner, seeds, seconds, nhits, nlass) and prints:

  1. per-rank work and time, with max/mean imbalance;
  2. how well seed count predicts task time (least-squares fit, R^2,
     Spearman rank correlation), per task and per rank;
  3. what a scheduler could gain: greedy longest-first (LPT) assignment of the
     measured tasks to the same number of ranks, using (a) the true task times
     (an oracle bound) and (b) seed counts as the cost estimate. Done at two
     granularities: unordered contig pairs (today's unit: (i,j), (j,i), N and C
     together) and single tasks (type, anchor, partner).

Section 3 is a model, not a measurement: it assumes a task takes the same time
on any rank and ignores everything outside search_seeds.

Section 4 (only if the task files carry the diagnostic columns ndiag, maxdiag,
nlcov, ...): how well each per-task feature predicts task time, a log-space
least-squares model on features known before alignment, the most expensive
tasks with their features, and scheduling bounds that use the best
pre-alignment predictor and a dynamic queue (model, not measurement).

Exit status 1 if an invariant fails (task seeds != reported seeds per rank).
Standard library only. Usage: analyze_align_balance.py <run-dir> [--csv out.csv]
"""
import glob
import math
import os
import re
import sys


def parse_kv(line):
    return {k: v for k, v in re.findall(r"(\w+)=(\S+)", line)}


def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def imbalance(xs):
    m = mean(xs)
    return max(xs) / m if xs and m > 0 else float("nan")


def linfit(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan"), float("nan"), float("nan")
    mx, my = mean(xs), mean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return float("nan"), float("nan"), float("nan")
    a = sxy / sxx
    b = my - a * mx
    r2 = (sxy * sxy) / (sxx * syy)
    return a, b, r2


def ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2.0
        i = j + 1
    return r


def spearman(xs, ys):
    _, _, r2 = linfit(ranks(xs), ranks(ys))
    return math.sqrt(r2) if r2 == r2 else float("nan")


def lpt(costs, true_times, nbins):
    """Assign tasks largest-cost first to the least-loaded bin (by cost);
    return the per-bin sums of true time."""
    order = sorted(range(len(costs)), key=lambda i: (-costs[i], i))
    load = [0.0] * nbins
    real = [0.0] * nbins
    for i in order:
        b = min(range(nbins), key=lambda k: (load[k], k))
        load[b] += costs[i]
        real[b] += true_times[i]
    return real


def lstsq(rows, ys):
    """Ordinary least squares via normal equations (small k). Returns coefs."""
    k = len(rows[0])
    A = [[sum(r[i] * r[j] for r in rows) for j in range(k)] for i in range(k)]
    for i in range(1, k):                   # tiny ridge on the slopes for stability
        A[i][i] *= 1.0 + 1e-9
    b = [sum(r[i] * y for r, y in zip(rows, ys)) for i in range(k)]
    for c in range(k):                      # Gauss-Jordan with partial pivoting
        p = max(range(c, k), key=lambda i: abs(A[i][c]))
        A[c], A[p] = A[p], A[c]
        b[c], b[p] = b[p], b[c]
        if abs(A[c][c]) < 1e-300:
            return [0.0] * k
        for i in range(k):
            if i != c:
                f = A[i][c] / A[c][c]
                A[i] = [a - f * x for a, x in zip(A[i], A[c])]
                b[i] -= f * b[c]
    return [b[i] / A[i][i] for i in range(k)]


def dynamic(order_costs, true_times, P):
    """List scheduling: tasks in the given order, each to the first free rank."""
    import heapq
    h = [0.0] * P
    for i in sorted(range(len(true_times)), key=lambda i: (-order_costs[i], i)):
        heapq.heapreplace(h, h[0] + true_times[i])
    return max(h)


def diagnose(tasks, feats, P):
    print("\n== 4. what predicts task time? (per task; model fits, not runs) ==")
    secs = [t[5] for t in tasks]
    num = lambda k: [float(v[k]) for v in feats]
    F = {"seeds": num("seeds"), "ndiag": num("ndiag"), "maxdiag": num("maxdiag"),
         "alen*blen": [float(v["alen"]) * float(v["blen"]) for v in feats],
         "nlass (after)": num("nlass"), "nlcov (after)": num("nlcov"),
         "nhits (after)": num("nhits")}
    keep = [i for i in range(len(secs)) if secs[i] > 0 and F["seeds"][i] > 0]
    ly = [math.log(secs[i]) for i in keep]
    print("%-15s %8s %9s %10s" % ("feature", "R^2", "spearman", "R^2(log)"))
    for name, xs in F.items():
        _, _, r2 = linfit(xs, secs)
        sp = spearman(xs, secs)
        lx = [math.log(xs[i] + 1.0) for i in keep]
        _, _, r2l = linfit(lx, ly)
        print("%-15s %8.3f %9.3f %10.3f" % (name, r2, sp, r2l))

    # log-space model on pre-alignment features only
    pre = []
    for k in ["seeds", "ndiag", "maxdiag", "alen*blen"]:
        lx = [math.log(F[k][i] + 1.0) for i in keep]
        if max(lx) - min(lx) > 1e-9:          # constant features make the fit singular
            pre.append(k)
    rows = [[1.0] + [math.log(F[k][i] + 1.0) for k in pre] for i in keep]
    coef = lstsq(rows, ly)
    pred_l = [sum(c * x for c, x in zip(coef, r)) for r in rows]
    myl = mean(ly)
    ss_res = sum((y - p) ** 2 for y, p in zip(ly, pred_l))
    ss_tot = sum((y - myl) ** 2 for y in ly)
    print("log-space model on pre-alignment features %s: R^2(log)=%.3f  coefs=%s"
          % (pre, 1 - ss_res / ss_tot if ss_tot else float("nan"),
             " ".join("%.3f" % c for c in coef)))
    pred = [0.0] * len(secs)
    for i, pl in zip(keep, pred_l):
        pred[i] = math.exp(pl)
    _, _, r2p = linfit(pred, secs)
    print("  -> predicted vs actual seconds: R^2=%.3f spearman=%.3f"
          % (r2p, spearman(pred, secs)))

    tot = sum(secs) or 1.0
    order = sorted(range(len(secs)), key=lambda i: -secs[i])
    print("most expensive tasks:")
    print("  %5s %4s %6s %6s %10s %8s %8s %10s %7s %12s %10s %10s"
          % ("share", "type", "anchor", "partn", "seeds", "ndiag", "maxdiag",
             "secs", "nlass", "nlcov", "alen", "blen"))
    for i in order[:12]:
        v = feats[i]
        print("  %4.1f%% %4s %6s %6s %10s %8s %8s %10.2f %7s %12s %10s %10s"
              % (100 * secs[i] / tot, v["type"], v["anchor"], v["partner"], v["seeds"],
                 v["ndiag"], v["maxdiag"], secs[i], v["nlass"], v["nlcov"],
                 v["alen"], v["blen"]))

    # scheduling bounds at the single-task level (model)
    m = tot / P
    est_seeds = F["seeds"]
    print("scheduling of single tasks on P=%d (slowest rank, x mean %.2f s):" % (P, m))
    for name, cost in (("true time (oracle)", secs), ("seeds", est_seeds),
                       ("pre-alignment model", pred)):
        st = max(lpt(cost, secs, P))
        dy = dynamic(cost, secs, P)
        print("  %-20s static LPT %.2f s (%.3fx)   dynamic queue %.2f s (%.3fx)"
              % (name, st, st / m, dy, dy / m))
    print("  largest single task %.2f s (floor for any scheduler)" % max(secs))


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    d = sys.argv[1]
    csv_out = sys.argv[3] if len(sys.argv) > 3 and sys.argv[2] == "--csv" else None

    ranks_kv = {}
    with open(os.path.join(d, "merge.err"), errors="replace") as f:
        for line in f:
            m = re.match(r"\[pr1\] rank (\d+): (.*)", line)
            if m:
                ranks_kv[int(m.group(1))] = parse_kv(m.group(2))
    if not ranks_kv:
        sys.exit("no [pr1] lines in %s/merge.err" % d)
    P = len(ranks_kv)

    tasks = []          # (rank, type, anchor, partner, seeds, secs, nlass)
    feats = []          # dict of all columns per task, aligned with tasks
    for path in sorted(glob.glob(os.path.join(d, "pairs.rank*.tsv"))):
        r = int(re.search(r"pairs\.rank(\d+)\.tsv", path).group(1))
        with open(path) as f:
            cols = next(f).rstrip("\n").split("\t")
            for line in f:
                v = dict(zip(cols, line.rstrip("\n").split("\t")))
                tasks.append((r, v["type"], int(v["anchor"]), int(v["partner"]),
                              int(v["seeds"]), float(v["seconds"]), int(v["nlass"])))
                feats.append(v)

    ok = True
    print("run: %s   ranks: %d   tasks: %d" % (d, P, len(tasks)))
    print("\n== 1. per rank ==")
    hdr = ("rank", "seeds", "task_s", "search_s", "pair_srch_s", "merge_core_s",
           "convert_s", "merge+align_s")
    print("%4s %12s %9s %9s %11s %12s %9s %13s" % hdr)
    cols = {k: [] for k in hdr[1:]}
    for r in sorted(ranks_kv):
        kv = ranks_kv[r]
        seeds = int(kv["seeds_n"]) + int(kv["seeds_c"])
        tsk = [t for t in tasks if t[0] == r]
        task_seeds = sum(t[4] for t in tsk)
        task_s = sum(t[5] for t in tsk)
        if task_seeds != seeds:
            print("  INVARIANT FAILED rank %d: task seeds %d != seeds_n+seeds_c %d"
                  % (r, task_seeds, seeds))
            ok = False
        search = float(kv["search_n_s"]) + float(kv["search_c_s"])
        pair = float(kv["pair_search_s"])
        merge_core = float(kv["merge_core_s"])
        conv = float(kv["convert_s"])
        tot = (float(kv["adapter_s"]) + merge_core + float(kv["redist_n_s"])
               + float(kv["redist_c_s"]) + float(kv["recompute_s"]) + pair + conv)
        row = (seeds, task_s, search, pair, merge_core, conv, tot)
        for k, v in zip(hdr[1:], row):
            cols[k].append(v)
        print("%4d %12d %9.3f %9.3f %11.3f %12.3f %9.3f %13.3f" % ((r,) + row))
    print("max/mean " + " ".join(
        "%s=%.3f" % (k, imbalance(cols[k])) for k in hdr[1:]))

    print("\n== 2. does seed count predict time? ==")
    xs = [t[4] for t in tasks]
    ys = [t[5] for t in tasks]
    a, b, r2 = linfit(xs, ys)
    print("per task (n=%d): secs = %.3e*seeds + %.3e   R^2=%.3f   spearman=%.3f"
          % (len(tasks), a, b, r2, spearman(xs, ys)))
    a2, b2, r2l = linfit(xs, [t[6] for t in tasks])
    print("per task: alignments(nlass) vs seeds R^2=%.3f" % r2l)
    big = sorted(tasks, key=lambda t: -t[5])[:5]
    tot_t = sum(ys) or 1.0
    print("largest tasks by time (share of all task time):")
    for t in big:
        print("  rank %d %s anchor %d partner %d: seeds=%d secs=%.3f (%.1f%%)"
              % (t[0], t[1], t[2], t[3], t[4], t[5], 100 * t[5] / tot_t))
    _, _, r2r = linfit(cols["seeds"], cols["search_s"])
    print("per rank (n=%d): search_s vs seeds R^2=%.3f" % (P, r2r))

    print("\n== 3. model: best-case scheduling of the measured tasks (prediction) ==")
    actual = max(cols["task_s"]) if cols["task_s"] else float("nan")
    lower = max(sum(ys) / P, max(ys) if ys else 0)
    print("actual max per-rank task time: %.3f s; lower bound max(sum/P, largest task): %.3f s"
          % (actual, lower))
    groups = {}
    for t in tasks:
        key = (min(t[2], t[3]), max(t[2], t[3]))
        g = groups.setdefault(key, [0, 0.0])
        g[0] += t[4]
        g[1] += t[5]
    for name, units in (("unordered pair", list(groups.values())),
                        ("single task", [[t[4], t[5]] for t in tasks])):
        seeds = [u[0] for u in units]
        times = [u[1] for u in units]
        oracle = max(lpt(times, times, P))
        byseed = max(lpt([float(s) for s in seeds], times, P))
        print("  %-15s units=%6d  LPT(true time) max=%.3f s (%.3f x mean)   "
              "LPT(seeds) max=%.3f s (%.3f x mean)"
              % (name, len(units), oracle, oracle / (sum(times) / P),
                 byseed, byseed / (sum(times) / P)))

    if feats and "maxdiag" in feats[0]:
        diagnose(tasks, feats, P)

    if csv_out:
        with open(csv_out, "w") as f:
            f.write("rank,type,anchor,partner,seeds,seconds,nlass\n")
            for t in tasks:
                f.write("%d,%s,%d,%d,%d,%.6f,%d\n" % t)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
