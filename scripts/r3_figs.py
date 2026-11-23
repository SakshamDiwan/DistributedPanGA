#!/usr/bin/env python3
"""R3 figures (PDF) from pga-prefix-hist output and existing run logs.

  r3_bucket_hist.pdf     human: records per 5-base bucket, ideal per-rank load at 64/128 ranks
  r3_rank_loads.pdf      human, 64 ranks: per-rank load / mean, production splitter vs optimal contiguous
  r3_balance_vs_ranks.pdf  max/mean rank load vs rank count: 5-base greedy (production), 5-base optimal,
                         8-base optimal; and the whole-contig extraction split
  r3_stage2_breakdown.pdf  critical-path decomposition of the sort stage (existing sweep logs)
  r3_sort_before_after.pdf sort stage and stage-2 memory, pre-R3 vs R3 defaults (bench round 3)

Usage: r3_figs.py --hist DIR [--results DIR] [--out report/fig]
Colours: reference categorical order (blue, orange, aqua, yellow, magenta), fixed by role.
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r3_bucket_analysis as BA          # noqa: E402
import r3_stage2_breakdown as SB         # noqa: E402

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, MUTED = "#1f1f1e", "#6f6e69"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.5, "legend.frameon": False,
                     "pdf.fonttype": 42})

def fig_hist(h10, out):
    total = h10.sum()
    fig, ax = plt.subplots(figsize=(6.5, 2.2))
    ax.bar(np.arange(1024), h10 / 1e6, width=1.0, color=C[0], linewidth=0)
    for P, ls in ((64, "--"), (128, ":")):
        ax.axhline(total / P / 1e6, color=INK, lw=0.8, ls=ls)
        ax.text(1030, total / P / 1e6, f"ideal load, {P} ranks", va="center", ha="left", fontsize=7, color=INK)
    for b in np.argsort(h10)[::-1][:2]:
        ax.annotate(BA.decode(int(b), 5), (b, h10[b] / 1e6), xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=7, color=INK)
    ax.set_xlim(-5, 1024 + 5)
    ax.set_xlabel("5-base prefix bucket (lexicographic, AAAAA = 0 ... TTTTT = 1023)")
    ax.set_ylabel("records (millions)")
    ax.set_title("hsapiens90.chr6: records per bucket", loc="left", fontsize=9, color=INK)
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)

def fig_rank_loads(h10, P, out):
    L = BA.loads_from_select(h10, BA.destination_map(h10, P), P)
    O = BA.optimal_contiguous(h10, P)
    fig, ax = plt.subplots(figsize=(6.5, 2.2))
    x = np.arange(P)
    ax.plot(x, L / L.mean(), color=C[0], lw=1.5, marker="o", ms=3, label=f"production splitter (max/mean {L.max() / L.mean():.3f})")
    ax.plot(x, O / O.mean(), color=C[1], lw=1.5, marker="s", ms=3, ls="--", label=f"optimal contiguous split (max/mean {O.max() / O.mean():.3f})")
    ax.axhline(1, color=MUTED, lw=0.6)
    ax.set_xlabel("rank"); ax.set_ylabel("records / mean")
    ax.set_title(f"hsapiens90.chr6, {P} ranks: sorted-slice size per rank (same 1024 buckets)", loc="left", fontsize=9, color=INK)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=2, fontsize=7)
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)

def fig_balance(h24, ctg, out):
    h10 = BA.coarsen(h24, 10); h16 = BA.coarsen(h24, 16)
    Ps = [8, 16, 32, 64, 96, 128, 256, 512]
    TICKS = [8, 16, 32, 64, 128, 256, 512]
    g = [BA.lstats(BA.loads_from_select(h10, BA.destination_map(h10, P), P))["max_over_mean"] for P in Ps]
    o5 = [BA.lstats(BA.optimal_contiguous(h10, P))["max_over_mean"] for P in Ps]
    o8 = [BA.lstats(BA.optimal_contiguous(h16, P))["max_over_mean"] for P in Ps]
    ex = []
    for P in Ps:
        _, R = BA.extraction_split(ctg, sorted(ctg), P); ex.append(R.max() / R.mean())
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.4))
    ax = axs[0]
    for y, c, m, ls, lab in ((g, C[0], "o", "-", "5-base, production"), (o5, C[1], "s", "--", "5-base, optimal contiguous"),
                             (o8, C[2], "^", ":", "8-base, optimal contiguous")):
        ax.plot(Ps, y, color=c, marker=m, ms=3.5, lw=1.5, ls=ls, label=lab)
    ax.set_xscale("log", base=2); ax.set_xticks(TICKS); ax.set_xticklabels(TICKS); ax.minorticks_off()
    ax.set_xlabel("ranks"); ax.set_ylabel("max / mean records per rank")
    ax.set_title("sort: slice balance", loc="left", fontsize=9, color=INK); ax.legend(fontsize=7)
    ax = axs[1]
    ax.plot(Ps, ex, color=C[0], marker="o", ms=3.5, lw=1.5)
    ax.set_xscale("log", base=2); ax.set_xticks(TICKS); ax.set_xticklabels(TICKS); ax.minorticks_off(); ax.set_yscale("log"); ax.set_yticks([1, 1.5, 2, 3, 4, 6]); ax.set_yticklabels(["1", "1.5", "2", "3", "4", "6"]); ax.minorticks_off()
    ax.set_xlabel("ranks"); ax.set_ylabel("max / mean records per rank")
    ax.set_title("extraction: whole-contig split", loc="left", fontsize=9, color=INK)
    fig.suptitle("hsapiens90.chr6", x=0.01, ha="left", fontsize=9, color=INK)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)

def fig_breakdown(results, out):
    import glob, statistics as st, re
    cfgs = [("hsapiens90", 4, 32), ("hsapiens90", 4, 64), ("primates14", 1, 16), ("primates14", 2, 32), ("primates14", 4, 64)]
    comps = [("hist_work_max", "histogram + Allreduce"), ("pack_s_max", "pack"), ("xfer_s", "exchange"),
             ("local_sort_s_max", "local sort + LCP"), ("verify_crit_s", "verify")]
    rows = []
    for ds, n, r in cfgs:
        vals = []
        for p in sorted(glob.glob(os.path.join(results, ds, f"n{n}-r{r}", "rep*", "merge.err"))):
            o1, mem, pr1, rs = SB.parse_run(p); s, _ = SB.run_stats(o1, mem, pr1, rs); vals.append(s)
        rows.append((f"{ds.replace('hsapiens90', 'human').replace('primates14', 'primates14')}\n{r} ranks / {n} nodes",
                     [st.median([v[k] for v in vals]) for k, _ in comps], st.median([v["sort_wall_s"] for v in vals])))
    fig, ax = plt.subplots(figsize=(6.5, 2.6))
    y = np.arange(len(rows))[::-1]
    left = np.zeros(len(rows))
    for i, (k, lab) in enumerate(comps):
        w = np.array([r[1][i] for r in rows])
        ax.barh(y, w, left=left, height=0.6, color=C[i], edgecolor="white", linewidth=1, label=lab)
        left += w
    for yi, r in zip(y, rows):
        ax.text(r[2] + 0.2, yi, f"{r[2]:.1f} s", va="center", fontsize=7, color=INK)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=7)
    ax.set_xlabel("seconds (slowest rank per substage; label = measured stage wall time)")
    ax.legend(ncol=5, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.45, 1.16))
    ax.grid(axis="y", visible=False)
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)

def fig_before_after(tsv, out):
    """Sort-stage wall time and per-rank stage-2 peak, pre-R3 vs final defaults (bench2 summary)."""
    import csv
    rows = {(r["dataset"], int(r["nodes"]), int(r["ranks"]), r["variant"]): r
            for r in csv.DictReader(open(tsv), delimiter="\t")}
    cfgs = [("hsapiens90", 4, 32), ("hsapiens90", 4, 64), ("hsapiens90", 8, 128), ("hsapiens90", 8, 256),
            ("primates14", 1, 16), ("primates14", 4, 64)]
    lab = [f"{'human' if d == 'hsapiens90' else d} {r}r/{n}n" for d, n, r in cfgs]
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.6))
    x = np.arange(len(cfgs)); w = 0.38
    for ax, key, ylab in ((axs[0], "sort_wall", "sort stage (s)"), (axs[1], "s2_peak", "stage-2 peak per rank (GB)")):
        pre = [float(rows[c + ("pre",)][key]) for c in cfgs]
        new = [float(rows[c + ("default",)][key]) for c in cfgs]
        ax.bar(x - w / 2, pre, w, color=C[1], label="pre-R3", edgecolor="white", linewidth=1)
        ax.bar(x + w / 2, new, w, color=C[0], label="R3 defaults", edgecolor="white", linewidth=1)
        ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=6.5, rotation=30, ha="right"); ax.set_ylabel(ylab)
        ax.grid(axis="x", visible=False)
    axs[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hist", required=True)
    ap.add_argument("--results", default=os.path.expanduser("~/pga-report-results-order"))
    ap.add_argument("--out", default="report/fig")
    ap.add_argument("--bench2", default=os.path.expanduser("~/pga-r3-results/bench2/summary.tsv"))
    a = ap.parse_args()
    h24, ctg = BA.load(a.hist, "hs")
    h10 = BA.coarsen(h24, 10)
    fig_hist(h10, os.path.join(a.out, "r3_bucket_hist.pdf"))
    fig_rank_loads(h10, 64, os.path.join(a.out, "r3_rank_loads.pdf"))
    fig_balance(h24, ctg, os.path.join(a.out, "r3_balance_vs_ranks.pdf"))
    fig_breakdown(a.results, os.path.join(a.out, "r3_stage2_breakdown.pdf"))
    if os.path.exists(a.bench2):
        fig_before_after(a.bench2, os.path.join(a.out, "r3_sort_before_after.pdf"))
    print("wrote", ", ".join(f"r3_{x}.pdf" for x in ("bucket_hist", "rank_loads", "balance_vs_ranks", "stage2_breakdown")))

if __name__ == "__main__":
    main()
