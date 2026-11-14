"""Human chr6 figure: time per stage, FastGA vs DistributedPanGA.

Numbers are copied from the run logs (hsc53, 2026-11-13/14); rerun this
script after editing them: `module load python; python3 make_human_fig.py`.

DistributedPanGA (4 nodes, FastGA 178afcb code, binary sha 4b96e0a85855f5e4):
  32 ranks: job 61255434 (repeat 61277457 within 1%)
  64 ranks: job 61244457 (repeat 61277456 within 1%)
  Extraction = slowest rank's extract_s; sort = max(extract_s + sort_stage_s)
  - extraction, so the wait for the slowest extractor counts as extraction
  (same rule as Saksham's tables). Seeds and alignment from the stage CSV.
FastGA (1 node): GIXmake -T16, job 61244453 (-T32 crashes on this input);
  FastGA + ALNtoPAF -T16 job 61255430, -T32 job 61244455.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = [  # label, [(stage, seconds), ...]
    ("FastGA, 16 threads", [("Index (GIXmake)", 400.39), ("Align (FastGA)", 1069.81), ("PAF", 10.22)]),
    ("FastGA, 32 threads", [("Index (GIXmake)", 400.39), ("Align (FastGA)", 1330.71), ("PAF", 5.83)]),
    ("DistributedPanGA, 32 ranks", [("Extraction", 9.86), ("Sort", 22.65), ("Seeds", 6.46), ("Alignment", 252.20)]),
    ("DistributedPanGA, 64 ranks", [("Extraction", 6.13), ("Sort", 18.16), ("Seeds", 3.85), ("Alignment", 151.92)]),
]
colors = {
    "Index (GIXmake)": "#9db4c0", "Align (FastGA)": "#5c6b73", "PAF": "#253237",
    "Extraction": "#f2c14e", "Sort": "#f78154", "Seeds": "#4d9078", "Alignment": "#b4436c",
}

fig, ax = plt.subplots(figsize=(7.0, 2.4))
seen = set()
for y, (label, parts) in enumerate(reversed(rows)):
    left = 0.0
    for stage, sec in parts:
        ax.barh(y, sec, left=left, color=colors[stage], edgecolor="white", linewidth=0.5,
                label=stage if stage not in seen else None)
        seen.add(stage)
        left += sec
    ax.text(left + 15, y, f"{left:,.0f} s", va="center", fontsize=8)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in reversed(rows)], fontsize=8)
ax.set_xlabel("Time (s)", fontsize=8)
ax.tick_params(axis="x", labelsize=8)
ax.set_xlim(0, 2000)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(fontsize=7, ncol=4, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.32))
fig.tight_layout()
fig.savefig("human_chr6.pdf", bbox_inches="tight")
