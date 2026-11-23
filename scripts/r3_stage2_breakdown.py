#!/usr/bin/env python3
"""R3: extraction / distributed-sort breakdown from existing run logs.

Reads the per-rank [o1-0] and [pr1-mem] lines in <results>/<ds>/n<N>-r<R>/rep<k>/merge.err
(layout of jobs/human/07_report_sweep.sbatch, PGA_RESULTS) and prints, per configuration
(median over repeats of each statistic):

  * extraction per rank: time, contigs, bases, records, idle ranks, max/mean;
  * stage-2 substages per rank: mean / max, and how much of each is waiting;
  * the sort-stage wall time with the extraction wait removed;
  * slice (post-exchange) records per rank: max/mean, CV, bytes moved;
  * per-stage memory (VmHWM per stage, reset at each stage boundary).

Waiting model (no barrier between extraction and sort in these builds):
  rank r enters run_stage2 at e_r and leaves the hist Allreduce no earlier than
  max_r(e_r). So hist_s(r) = count_r + (max_e - e_r) + allreduce + map, i.e.
  hist work ~= hist_s(r) - (max_e - e_r). local_sort_s and pack_s contain no
  collective, so their per-rank values are pure computation; exchange_s and
  verify_s start with a collective and absorb the spread of the preceding
  local phase.

Usage: r3_stage2_breakdown.py [RESULTS_DIR] [--tsv OUT.tsv] [--ranks OUTDIR]
"""
import os, re, sys, glob, statistics as st

def kv(line):
    return {a: b for a, b in (x.split("=", 1) for x in line.split() if "=" in x)}

def parse_run(path):
    o1, mem, pr1 = {}, {}, {}
    rsize = None
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("[o1-0] rank"):
                r = int(line.split()[2].rstrip(":")); o1[r] = {k: float(v) for k, v in kv(line).items()}
            elif line.startswith("[pr1-mem] rank"):
                r = int(line.split()[2].rstrip(":")); mem[r] = {k: float(v) for k, v in kv(line).items()}
            elif line.startswith("[pr1] rank"):
                r = int(line.split()[2].rstrip(":")); d = kv(line)
                pr1[r] = d
                if rsize is None and float(d.get("slice_records", 0)) > 0:
                    rsize = float(d["record_bytes"]) / float(d["slice_records"])
    return o1, mem, pr1, rsize

def cv(xs):
    m = st.mean(xs)
    return st.pstdev(xs) / m if m else float("nan")

def run_stats(o1, mem, pr1, rsize):
    R = sorted(o1)
    n = len(R)
    g = lambda k: [o1[r][k] for r in R]
    e = g("extract_s"); emax = max(e)
    s = {}
    s["ranks"] = n
    s["rsize"] = rsize or float("nan")
    rec = g("records"); sl = g("slice_records")
    tot = sum(rec)
    s["records_total"] = tot
    assert abs(sum(sl) - tot) < 0.5, "slice and extracted totals differ"
    # extraction
    s["ext_max_s"] = emax; s["ext_mean_s"] = st.mean(e); s["ext_min_s"] = min(e)
    s["ext_idle_ranks"] = sum(1 for x in rec if x == 0)
    s["ext_rec_max_over_mean"] = max(rec) / st.mean(rec)
    s["ext_bases_max_over_mean"] = max(g("bases")) / st.mean(g("bases"))
    s["ext_time_max_over_mean"] = emax / st.mean(e)
    s["ext_wait_mean_s"] = st.mean(emax - x for x in e)
    eb = [b for b, t in zip(g("bases"), e) if t > 0 and b > 0]
    s["ext_Mbp_per_s_median"] = st.median([b / t / 1e6 for b, t in zip(g("bases"), e) if t > 0 and b > 0]) if eb else float("nan")
    s["fallback_reads_total"] = sum(g("pack_reads"))
    s["window_packs_frac"] = sum(g("window_packs")) / tot if tot else float("nan")
    # sort stage
    ss = g("sort_stage_s")
    s["sort_raw_max_s"] = max(ss)
    s["sort_wall_s"] = max(a + b for a, b in zip(e, ss)) - emax  # from last extractor to last sorter
    for k in ("hist_s", "pack_s", "exchange_s", "local_sort_s", "verify_s"):
        v = g(k); s[k + "_mean"] = st.mean(v); s[k + "_max"] = max(v)
    hw = [h - (emax - x) for h, x in zip(g("hist_s"), e)]
    s["hist_work_max"] = max(hw); s["hist_work_mean"] = st.mean(hw)
    # critical-path decomposition of sort_wall:
    #   hist (after last extractor) + pack (slowest) + exchange (after slowest packer)
    #   + local sort (slowest) + verify (after slowest sorter)
    # exchange on the slowest-packer rank ~= transfer; verify on the slowest sorter ~= true verify.
    ip = max(range(n), key=lambda i: g("pack_s")[i])
    il = max(range(n), key=lambda i: g("local_sort_s")[i])
    s["xfer_s"] = g("exchange_s")[ip]               # exchange time seen by the slowest packer
    s["verify_crit_s"] = g("verify_s")[il]          # verify time seen by the slowest sorter
    s["sync_wait_frac"] = 1 - (s["hist_work_max"] + s["pack_s_max"] + s["xfer_s"]
                               + s["local_sort_s_max"] + s["verify_crit_s"]) / s["sort_wall_s"] \
        if s["sort_wall_s"] > 0 else float("nan")
    s["comm_frac"] = s["xfer_s"] / s["sort_wall_s"] if s["sort_wall_s"] > 0 else float("nan")
    s["localsort_frac"] = s["local_sort_s_max"] / s["sort_wall_s"] if s["sort_wall_s"] > 0 else float("nan")
    # slice balance
    s["slice_max"] = max(sl); s["slice_mean"] = st.mean(sl); s["slice_min"] = min(sl)
    s["slice_max_over_mean"] = max(sl) / st.mean(sl); s["slice_cv"] = cv(sl)
    s["slice_empty"] = sum(1 for x in sl if x == 0)
    s["slice_under_half"] = sum(1 for x in sl if x < 0.5 * st.mean(sl))
    s["sort_Mrec_per_s_median"] = st.median([a / b / 1e6 for a, b in zip(sl, g("local_sort_s")) if b > 0])
    if rsize:
        s["bytes_moved_GB"] = tot * rsize / 1e9      # every record is sent (incl. to self)
        s["recv_max_GB"] = max(sl) * rsize / 1e9
    # memory (kB -> GB)
    if mem:
        M = sorted(mem)
        for st_ in ("extract", "stage2", "adapter", "merge", "redist", "align", "convert"):
            v = [mem[r][st_ + "_peak_kb"] / 1e6 for r in M]
            s["mem_" + st_ + "_max_GB"] = max(v); s["mem_" + st_ + "_mean_GB"] = st.mean(v)
        ov = [mem[r]["overall_peak_kb"] / 1e6 for r in M]
        s["mem_overall_max_GB"] = max(ov); s["mem_overall_sum_GB"] = sum(ov)
        s["mem_reset_ok"] = min(mem[r]["reset_ok"] for r in M)
        s["mem_peak_stage2_frac_ranks"] = sum(1 for r in M if mem[r]["stage2_peak_kb"] == mem[r]["overall_peak_kb"]) / len(M)
    return s, {r: dict(o1[r], **{"mem_" + k: v for k, v in mem.get(r, {}).items()}) for r in R}

def main():
    a = sys.argv[1:]
    root = a[0] if a and not a[0].startswith("--") else os.path.expanduser("~/pga-report-results-order")
    tsv = a[a.index("--tsv") + 1] if "--tsv" in a else None
    rdir = a[a.index("--ranks") + 1] if "--ranks" in a else None
    cfg = {}
    for p in sorted(glob.glob(os.path.join(root, "*", "n*-r*", "rep*", "merge.err"))):
        parts = p.split(os.sep)
        ds, nr, rep = parts[-4], parts[-3], parts[-2]
        m = re.match(r"n(\d+)-r(\d+)", nr)
        o1, mem, pr1, rsize = parse_run(p)
        if not o1:
            continue
        s, per = run_stats(o1, mem, pr1, rsize)
        cfg.setdefault((ds, int(m.group(1)), int(m.group(2))), []).append((rep, s, per))
    keys = sorted(cfg, key=lambda k: (k[0], k[2], k[1]))
    cols = None
    rows = []
    for k in keys:
        reps = cfg[k]
        cols = cols or list(reps[0][1].keys())
        med = {c: st.median([r[1][c] for r in reps]) for c in cols}
        spread = {c: (max(r[1][c] for r in reps) - min(r[1][c] for r in reps)) for c in cols}
        rows.append((k, len(reps), med, spread))
        if rdir:
            os.makedirs(rdir, exist_ok=True)
            rep, _, per = sorted(reps, key=lambda r: r[0])[0]
            with open(os.path.join(rdir, f"{k[0]}-n{k[1]}-r{k[2]}-{rep}.ranks.tsv"), "w") as f:
                ks = ["records", "slice_records", "bases", "contigs", "extract_s", "hist_s", "pack_s",
                      "exchange_s", "local_sort_s", "verify_s", "sort_stage_s",
                      "mem_extract_peak_kb", "mem_stage2_peak_kb", "mem_adapter_peak_kb", "mem_overall_peak_kb"]
                f.write("rank\t" + "\t".join(ks) + "\n")
                for r in sorted(per):
                    f.write(f"{r}\t" + "\t".join(f"{per[r].get(x, float('nan')):g}" for x in ks) + "\n")
    if tsv:
        with open(tsv, "w") as f:
            f.write("dataset\tnodes\tranks\treps\t" + "\t".join(cols) + "\t" + "\t".join("spread_" + c for c in cols) + "\n")
            for k, n, med, spr in rows:
                f.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{n}\t" + "\t".join(f"{med[c]:.6g}" for c in cols)
                        + "\t" + "\t".join(f"{spr[c]:.4g}" for c in cols) + "\n")
    # human-readable summary
    show = ["records_total", "ext_max_s", "ext_idle_ranks", "ext_rec_max_over_mean", "ext_wait_mean_s",
            "sort_raw_max_s", "sort_wall_s", "hist_work_max", "pack_s_max", "xfer_s", "local_sort_s_max",
            "verify_crit_s", "sync_wait_frac", "comm_frac", "localsort_frac", "slice_max_over_mean", "slice_cv",
            "bytes_moved_GB", "recv_max_GB", "sort_Mrec_per_s_median", "mem_extract_max_GB", "mem_stage2_max_GB",
            "mem_adapter_max_GB", "mem_overall_max_GB", "mem_overall_sum_GB"]
    w = 12
    print("config".ljust(24) + "".join(c[:w - 1].rjust(w) for c in show))
    for k, n, med, spr in rows:
        print(f"{k[0]} n{k[1]} r{k[2]} x{n}".ljust(24) + "".join(f"{med.get(c, float('nan')):{w}.4g}" for c in show))

if __name__ == "__main__":
    main()
