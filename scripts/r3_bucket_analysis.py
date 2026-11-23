#!/usr/bin/env python3
"""R3 P0-A / P1-A / P1-B / P1-C: prefix-bucket and extraction balance from pga-prefix-hist output.

Inputs (per dataset): <dir>/<tag>[.cK].hist24.bin (int64[2^24], summed over chunks) and
<dir>/<tag>[.cK].contigs.tsv. Optional: a results tree (07_report_sweep layout) to check the
simulated per-rank loads against the logged slice_records / records of real runs.

Prints, and writes <out>/<tag>.*.tsv for figures:
  * the 1024-bucket histogram and its statistics, largest buckets with their 5-base prefixes;
  * per-rank loads under the production splitter (exact port of compute_destination_map);
  * alternatives that keep whole 12-base groups together: optimal contiguous split of the
    1024 buckets, contiguous splits at 6/8/12-base prefixes, LPT over buckets (non-contiguous);
  * msd_sort thread balance per rank (list scheduling of parts over 4 threads, 256 vs 1024 parts);
  * extraction balance of the whole-contig work split vs the ideal.

Usage: r3_bucket_analysis.py --hist DIR --tag hs --name hsapiens90 [--results DIR] [--out DIR]
"""
import argparse, glob, os, statistics as st, sys
import numpy as np

BASES = "ACGT"

def decode(prefix, nbases):
    return "".join(BASES[(prefix >> (2 * (nbases - 1 - i))) & 3] for i in range(nbases))

def load(dirn, tag):
    files = sorted(glob.glob(os.path.join(dirn, tag + ".hist24.bin")) + glob.glob(os.path.join(dirn, tag + ".c*.hist24.bin")))
    if not files:
        sys.exit(f"no hist24 for {tag} in {dirn}")
    h = np.zeros(1 << 24, dtype=np.int64)
    ctg = {}
    for f in files:
        h += np.fromfile(f, dtype="<i8")
        for line in open(f.replace(".hist24.bin", ".contigs.tsv")):
            c, clen, rec = map(int, line.split())
            ctg[c] = (clen, ctg.get(c, (0, 0))[1] + rec)
    return h, ctg

def coarsen(h24, bits):
    return h24.reshape(1 << bits, -1).sum(axis=1)

# --- exact port of src/dsort.c compute_destination_map -------------------------------
def destination_map(counts, P):
    """Exact port of compute_destination_map, O(P log n): at most one split per bucket.

    Returns select as an int array (bucket -> rank)."""
    cum = np.cumsum(np.asarray(counts, dtype=np.int64))
    n_b = len(cum)
    total = int(cum[-1])
    select = np.zeros(n_b, dtype=np.int64)
    n, t, i = 1, total // P, 0
    while i < n_b:
        # next bucket at or after i whose cumulative count reaches t
        j = int(np.searchsorted(cum, t, side="left"))
        if j < i: j = i
        if j >= n_b:
            select[i:] = n - 1
            break
        select[i:j] = n - 1
        prev = int(cum[j - 1]) if j > 0 else 0
        select[j] = n if (int(cum[j]) - t > t - prev) else n - 1
        n += 1
        if n > P: n = P
        t = (n * total) // P
        i = j + 1
    return select

def loads_from_select(counts, select, P):
    return np.bincount(np.asarray(select), weights=np.asarray(counts, dtype=np.float64), minlength=P).astype(np.int64)

def optimal_contiguous(counts, P):
    """Minimise the max load over contiguous splits into <= P parts (binary search + greedy jumps)."""
    c = np.asarray(counts, dtype=np.int64)
    cum = np.cumsum(c)
    n_b = len(c)
    def cut(cap):
        ends, start = [], 0
        while start < n_b:
            base = int(cum[start - 1]) if start else 0
            e = int(np.searchsorted(cum, base + cap, side="right"))   # first index exceeding cap
            if e <= start: return None
            ends.append(e); start = e
            if len(ends) > P: return None
        return ends
    lo, hi = int(c.max()), int(cum[-1])
    while lo < hi:
        mid = (lo + hi) // 2
        if cut(mid) is not None: hi = mid
        else: lo = mid + 1
    ends = cut(lo)
    L, prev = [], 0
    for e in ends:
        L.append(int(cum[e - 1]) - prev); prev = int(cum[e - 1])
    L += [0] * (P - len(L))
    return np.array(L, dtype=np.int64)

def lpt(counts, P):
    import heapq
    hp = [(0, r) for r in range(P)]
    for x in sorted((int(v) for v in counts), reverse=True):
        l, r = heapq.heappop(hp); heapq.heappush(hp, (l + x, r))
    return np.array(sorted(l for l, _ in hp), dtype=np.int64)

def lstats(L):
    m = L.mean()
    return dict(max=int(L.max()), min=int(L.min()), mean=m, max_over_mean=L.max() / m,
                cv=L.std() / m, empty=int((L == 0).sum()), under_half=int((L < 0.5 * m).sum()))

def msd_schedule(part_sizes, threads):
    """msd_sort: parts launched in index order on the next free thread (list scheduling)."""
    import heapq
    hp = [0.0] * threads
    for x in part_sizes:
        if x > 0:
            t = heapq.heappop(hp); heapq.heappush(hp, t + x)
    total = sum(part_sizes)
    mk = max(hp) if total else 0
    return mk, (total / threads) / mk if mk else 1.0

def extraction_split(ctg, ncontig_order, P):
    """compute_work_split: whole contigs, thresholds on cumulative bases (port of work_split.c)."""
    clen = [ctg[c][0] for c in ncontig_order]
    rec = [ctg[c][1] for c in ncontig_order]
    total = sum(clen)
    split = [0] * (P + 1)
    so_far, nxt, w = 0, total // P, 1
    for r, L in enumerate(clen):
        so_far += L
        while so_far >= nxt and w < P:
            split[w] = r + 1; w += 1; nxt = (total * w) // P
    split[P] = len(clen)
    B = np.array([sum(clen[split[i]:split[i + 1]]) for i in range(P)], dtype=np.int64)
    R = np.array([sum(rec[split[i]:split[i + 1]]) for i in range(P)], dtype=np.int64)
    return B, R

def logged_runs(results, name):
    out = {}
    for p in glob.glob(os.path.join(results, name, "n*-r*", "rep1", "merge.err")):
        sl, ex = {}, {}
        for line in open(p, errors="replace"):
            if line.startswith("[o1-0] rank"):
                r = int(line.split()[2].rstrip(":"))
                kv = dict(x.split("=", 1) for x in line.split() if "=" in x)
                sl[r] = int(kv["slice_records"]); ex[r] = int(kv["records"])
        if sl:
            P = len(sl)
            out[P] = (np.array([sl[r] for r in range(P)]), np.array([ex[r] for r in range(P)]))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hist", required=True); ap.add_argument("--tag", required=True)
    ap.add_argument("--name", required=True); ap.add_argument("--results", default=os.path.expanduser("~/pga-report-results-order"))
    ap.add_argument("--out", default=None); ap.add_argument("--ranks", default="8,16,32,64,96,127,128,256,512")
    a = ap.parse_args()
    out = a.out or a.hist
    os.makedirs(out, exist_ok=True)
    h24, ctg = load(a.hist, a.tag)
    h10 = coarsen(h24, 10)
    total = int(h10.sum())
    print(f"== {a.name}: {total:,} records, {sum(v[0] for v in ctg.values()):,} bases, {len(ctg)} contigs")

    # ---- 1024-bucket histogram ----
    nz = h10[h10 > 0]
    print(f"buckets: nonempty {len(nz)}/1024  max {h10.max():,}  min {h10.min():,}  mean {h10.mean():,.0f}  "
          f"median {int(np.median(h10)):,}  max/mean {h10.max() / h10.mean():.2f}  cv {h10.std() / h10.mean():.3f}  "
          f"max/min {h10.max() / max(1, h10.min()):.1f}")
    top = np.argsort(h10)[::-1][:10]
    print("largest buckets: " + ", ".join(f"{decode(b, 5)}({b}) {h10[b] / total * 100:.3f}%" for b in top))
    bot = np.argsort(h10)[:5]
    print("smallest buckets: " + ", ".join(f"{decode(b, 5)}({b}) {h10[b] / total * 100:.4f}%" for b in bot))
    with open(os.path.join(out, f"{a.tag}.buckets.tsv"), "w") as f:
        f.write("bucket\tprefix\tcount\n")
        for b in range(1024):
            f.write(f"{b}\t{decode(b, 5)}\t{h10[b]}\n")
    # largest 12-base group (the indivisible unit for any group-preserving split)
    g = int(h24.max())
    print(f"largest 12-base group: {decode(int(h24.argmax()), 12)} {g:,} records ({g / total * 100:.4f}% of all)")

    # ---- validate the splitter port against real runs ----
    runs = logged_runs(a.results, a.name)
    for P, (sl, ex) in sorted(runs.items()):
        L = loads_from_select(h10, destination_map(h10, P), P)
        ok = np.array_equal(L, sl)
        print(f"check P={P}: simulated slice loads == logged slice_records: {ok}")
        if not ok:
            print("   max diff", np.abs(L - sl).max())

    # ---- per-rank loads: production vs alternatives ----
    rows = []
    order = sorted(ctg)
    hdr = ("P", "scheme", "max_over_mean", "cv", "empty", "under_half", "largest_bucket_over_ideal")
    print("\n" + "\t".join(hdr))
    for P in [int(x) for x in a.ranks.split(",")]:
        ideal = total / P
        schemes = []
        sel = destination_map(h10, P)
        L = loads_from_select(h10, sel, P)
        schemes.append(("prod-5base-greedy", L, h10.max()))
        with open(os.path.join(out, f"{a.tag}.P{P}.rankloads.tsv"), "w") as f:
            f.write("rank\trecords\tfirst_bucket\tlast_bucket\tnbuckets\n")
            for r in range(P):
                bs = list(np.nonzero(sel == r)[0])
                f.write(f"{r}\t{L[r]}\t{bs[0] if bs else -1}\t{bs[-1] if bs else -1}\t{len(bs)}\n")
        schemes.append(("5base-optimal-contig", optimal_contiguous(h10, P), h10.max()))
        schemes.append(("5base-LPT-noncontig", lpt(h10, P), h10.max()))
        for bits, nb in ((12, 6), (16, 8), (24, 12)):
            hc = coarsen(h24, bits)
            schemes.append((f"{nb}base-greedy", loads_from_select(hc, destination_map(hc, P), P), hc.max()))
            schemes.append((f"{nb}base-optimal-contig", optimal_contiguous(hc, P), hc.max()))
        for name, L, big in schemes:
            s = lstats(L)
            rows.append((P, name, s, big / ideal))
            print(f"{P}\t{name}\t{s['max_over_mean']:.4f}\t{s['cv']:.4f}\t{s['empty']}\t{s['under_half']}\t{big / ideal:.3f}")
    with open(os.path.join(out, f"{a.tag}.schemes.tsv"), "w") as f:
        f.write("P\tscheme\tmax_over_mean\tcv\tempty\tunder_half\tlargest_unit_over_ideal\n")
        for P, name, s, b in rows:
            f.write(f"{P}\t{name}\t{s['max_over_mean']:.5f}\t{s['cv']:.5f}\t{s['empty']}\t{s['under_half']}\t{b:.4f}\n")

    # ---- msd_sort thread balance (production splitter) ----
    print("\nmsd_sort list-scheduling efficiency (mean over ranks; worst rank), work ~ records per part")
    h8 = coarsen(h24, 8)
    with open(os.path.join(out, f"{a.tag}.msdsched.tsv"), "w") as f:
        f.write("P\tthreads\tparts\tmean_eff\tmin_eff\tmean_nonempty_parts\tslowest_rank_rel\n")
        for P in [int(x) for x in a.ranks.split(",")]:
            sel = destination_map(h10, P)
            for threads in (2, 4, 8, 16):
                for parts in (256, 1024):
                    effs, nps, mks = [], [], []
                    for r in range(P):
                        bs = list(np.nonzero(sel == r)[0])
                        if not bs:
                            continue
                        if parts == 1024:
                            ps = [int(h10[b]) for b in bs]
                        else:
                            fb = sorted(set(b >> 2 for b in bs))
                            ps = [int(sum(h10[b] for b in bs if b >> 2 == x)) for x in fb]
                        mk, e = msd_schedule(ps, threads)
                        effs.append(e); nps.append(sum(1 for x in ps if x)); mks.append(mk)
                    rel = max(mks) / (total / P / threads)
                    f.write(f"{P}\t{threads}\t{parts}\t{st.mean(effs):.4f}\t{min(effs):.4f}\t{st.mean(nps):.2f}\t{rel:.4f}\n")
                    if threads == 4:
                        print(f"P={P:4d} threads=4 parts={parts:4d}: mean eff {st.mean(effs):.3f}  worst {min(effs):.3f}  "
                              f"parts/rank {st.mean(nps):.1f}  slowest rank makespan / ideal {rel:.3f}")

    # ---- extraction balance ----
    print("\nextraction (whole-contig split by bases): records max/mean, idle ranks, max-contig/ideal")
    maxc = max(v[1] for v in ctg.values())
    with open(os.path.join(out, f"{a.tag}.extract.tsv"), "w") as f:
        f.write("P\trec_max_over_mean\tbases_max_over_mean\tidle\tlargest_contig_over_ideal\n")
        for P in [int(x) for x in a.ranks.split(",")]:
            B, R = extraction_split(ctg, order, P)
            idle = int((R == 0).sum())
            f.write(f"{P}\t{R.max() / R.mean():.4f}\t{B.max() / B.mean():.4f}\t{idle}\t{maxc / (total / P):.4f}\n")
            chk = ""
            if P in runs:
                chk = "  (logged records per rank equal: %s)" % np.array_equal(R, runs[P][1])
            print(f"P={P:4d}: records max/mean {R.max() / R.mean():.3f}  idle {idle:3d}  largest contig / ideal {maxc / (total / P):.2f}{chk}")

if __name__ == "__main__":
    main()
