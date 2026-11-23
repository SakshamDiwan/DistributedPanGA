#!/usr/bin/env python3
"""Summarise pga-dsort-bench runs (jobs/r3/r3_dsort_bench.sbatch output).

Reads <root>/<dataset>/n<N>-r<R>/<variant>.rep<k>.err and prints, per configuration and variant,
medians over repeats (min-max in brackets) of:
  sort_wall   stage wall time between the barriers around run_stage2 (s)
  pack, xfer, rebucket, msd, lcp   slowest rank's substage time (s), from [r3-s2]
  sizes       slowest packer wait in the size Alltoall is excluded; xfer is the transfer itself
  s2_peak     max over ranks of the stage-2 VmHWM (GB);  sum_peak: sum over ranks of max(extract, stage 2) (GB)
and the digest check (kseq / multi must have one value per dataset across all variants and rank counts).

Usage: r3_bench_summary.py [ROOT] [--tsv OUT]
"""
import glob, os, re, sys, statistics as st
from collections import defaultdict

def kv(line):
    return {a: b for a, b in (x.split("=", 1) for x in line.split() if "=" in x)}

def parse(path):
    ranks, bench = [], None
    for line in open(path, errors="replace"):
        if line.startswith("[r3-s2] rank"):
            ranks.append({k: float(v) for k, v in kv(line).items() if re.match(r"^-?[0-9.e+-]+$", v)})
        elif line.startswith("[bench] rep="):
            bench = kv(line)
    return ranks, bench

def main():
    a = sys.argv[1:]
    root = a[0] if a and not a[0].startswith("--") else os.path.expanduser("~/pga-r3-results")
    tsv = a[a.index("--tsv") + 1] if "--tsv" in a else None
    R = defaultdict(list)       # (ds, n, r, variant) -> list of dict
    digests = defaultdict(set)  # ds -> {(kseq, multi)}
    for p in glob.glob(os.path.join(root, "*", "n*-r*", "*.rep*.err")):
        ds, cfg = p.split(os.sep)[-3], p.split(os.sep)[-2]
        if not re.match(r"n\d+-r\d+$", cfg):
            continue
        n, r = map(int, re.findall(r"\d+", cfg))
        var = os.path.basename(p).split(".rep")[0]
        ranks, b = parse(p)
        if not b or not ranks:
            continue
        digests[ds].add((b["kseq"], b["multi"]))
        mx = lambda k: max(x.get(k, 0) for x in ranks)
        R[(ds, n, r, var)].append(dict(
            sort_wall=float(b["sort_wall_s"]), extract_wall=float(b["extract_wall_s"]),
            pack=mx("pack_s"), xfer=mx("xfer_s"), rebucket=mx("rebucket_s"), msd=mx("msd_s"), lcp=mx("lcp_s"),
            s2_peak=float(b["max_s2_hwm_kb"]) / 1e6, sum_peak=float(b["sum_peak_kb"]) / 1e6,
            rounds=mx("rounds"), staging=mx("staging") / 1e9))
    cols = ["sort_wall", "pack", "xfer", "rebucket", "msd", "lcp", "s2_peak", "sum_peak", "rounds", "staging"]
    order = ["old", "take", "p1024", "tphys", "p1024tphys", "large", "r256M", "r4G", "opt", "all"]
    out = []
    print(f"{'config':22s}{'variant':12s}{'n':>3s}" + "".join(f"{c:>18s}" for c in cols))
    for key in sorted({k[:3] for k in R}, key=lambda k: (k[0], k[2], k[1])):
        base = R.get(key + ("take",))
        bmed = st.median(x["sort_wall"] for x in base) if base else None
        for v in sorted({k[3] for k in R if k[:3] == key}, key=lambda v: order.index(v) if v in order else 99):
            rs = R[key + (v,)]
            cells = []
            row = dict(dataset=key[0], nodes=key[1], ranks=key[2], variant=v, reps=len(rs))
            for c in cols:
                xs = [x[c] for x in rs]
                m = st.median(xs)
                row[c] = m; row[c + "_min"] = min(xs); row[c + "_max"] = max(xs)
                cells.append(f"{m:8.3f}[{min(xs):.2f}-{max(xs):.2f}]" if c not in ("rounds",) else f"{m:18.0f}")
            if bmed:
                row["sort_vs_take"] = row["sort_wall"] / bmed
            out.append(row)
            print(f"{key[0][:10]} n{key[1]} r{key[2]:<5d}".ljust(22) + f"{v:12s}{len(rs):3d}" + "".join(f"{c:>18s}" for c in cells))
    print("\ndigests (kseq, multi) per dataset:")
    for ds, s in sorted(digests.items()):
        print(f"  {ds}: {len(s)} distinct {'OK' if len(s) == 1 else 'MISMATCH'} {sorted(s)}")
    if tsv:
        keys = list(out[0].keys())
        with open(tsv, "w") as f:
            f.write("\t".join(keys) + "\n")
            for r in out:
                f.write("\t".join(f"{r.get(k, '')}" if not isinstance(r.get(k), float) else f"{r[k]:.6g}" for k in keys) + "\n")

if __name__ == "__main__":
    main()
