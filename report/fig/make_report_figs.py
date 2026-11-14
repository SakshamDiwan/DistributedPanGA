"""Report figures: strong_scaling.pdf and vs_fastga.pdf, all four datasets.

Reads the numbers from the table files next to this script, so the figures
and the tables cannot disagree:
  table_pipeline.tex  dataset, nodes, ranks, extraction, sort, seeds, alignment, total
  table_fastga.tex    dataset, threads, GIXmake, FastGA, ALNtoPAF, total
  table_perrank.tex   dataset, nodes, ranks, slowest, average, fastest, slowest/avg

Run:  module load python; python3 make_report_figs.py

Colours, fonts and line styles match the figures first drawn for the yeast
and primate datasets.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE = os.path.dirname(os.path.abspath(__file__))
BLUE, ORANGE, GREEN, AMBER = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
GREY, DARK, GRID = "#6b6b6b", "#222222", "#e3e3e3"

DATASETS = ["yeast7", "primates4.chr6", "primates14.chr6", "hsapiens90.chr6"]
DATASET_STYLE = dict(zip(DATASETS, [(BLUE, "o"), (ORANGE, "s"), (GREEN, "^"), (AMBER, "D")]))
STAGES = [("extraction", BLUE, "o"), ("distributed sort", ORANGE, "s"),
          ("seed generation + redistribution", GREEN, "^"), ("alignment + PAF output", AMBER, "D")]

plt.rcParams.update({
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "text.color": DARK, "axes.labelcolor": DARK, "axes.titlecolor": DARK,
    "axes.edgecolor": GREY, "xtick.color": GREY, "ytick.color": GREY,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "lines.linewidth": 1.5, "lines.markersize": 5, "legend.frameon": False,
})


def read_table(name):
    """Rows of a generated table file as lists: dataset string, then floats."""
    rows = []
    with open(os.path.join(HERE, name)) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("%") or line.startswith("\\"):
                continue
            cells = [c.strip() for c in line.rstrip("\\").split("&")]
            rows.append([cells[0]] + [float(c.replace(",", "")) for c in cells[1:]])
    return rows


def by_dataset(rows):
    out = {d: [] for d in DATASETS}
    for r in rows:
        out[r[0]].append(r[1:])
    return out


def plain(v, _pos=None):
    return f"{v:g}"


def log_axes(ax, xs, ymin, ymax):
    """Log-log axes with plain-number ticks at 1-2-5 steps."""
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.xaxis.set_major_locator(FixedLocator(xs))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(plain))
    ax.set_xlim(xs[0] / 1.15, xs[-1] * 1.15)
    lo, hi = ymin / 1.25, ymax * 1.25
    ticks = [m * 10 ** e for e in range(-2, 5) for m in (1, 2, 5)]
    ticks = [t for t in ticks if lo <= t <= hi]
    if len(ticks) < 3:  # narrow range: fall back to finer steps
        ticks = [m * 10 ** e for e in range(-2, 5) for m in (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8)]
        ticks = [t for t in ticks if lo <= t <= hi]
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(plain))
    ax.set_ylim(lo, hi)


def node_labels(ax, runs):
    """Rank tick labels; add the node count when a dataset uses several nodes."""
    if any(r[0] > 1 for r in runs):
        ax.set_xticklabels([f"{int(r[1])}\n({int(r[0])} node{'s' if r[0] > 1 else ''})" for r in runs])


def strong_scaling():
    pipe = by_dataset(read_table("table_pipeline.tex"))   # nodes, ranks, 4 stages, total
    rank = by_dataset(read_table("table_perrank.tex"))    # nodes, ranks, slow, avg, fast, ratio
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.7))
    a, panels, f = axes[0][0], [axes[0][1], axes[0][2], axes[1][0], axes[1][1]], axes[1][2]

    # (a) total time per dataset, with ideal scaling from the smallest run
    for d in DATASETS:
        c, m = DATASET_STYLE[d]
        x = [r[1] for r in pipe[d]]
        y = [r[6] for r in pipe[d]]
        a.plot(x, y, color=c, marker=m, label=d)
        a.plot(x, [y[0] * x[0] / v for v in x], color=c, ls=":", lw=1.0)
    every = sorted({r[1] for d in DATASETS for r in pipe[d]})
    totals = [r[6] for d in DATASETS for r in pipe[d]]
    log_axes(a, every, min(totals), max(totals))
    a.set_title("(a) Strong scaling (dotted: ideal)", loc="left")
    a.set_xlabel("MPI ranks"); a.set_ylabel("total time (s)")

    # (b)-(e) time per stage for each dataset
    for ax, d, tag in zip(panels, DATASETS, "bcde"):
        x = [r[1] for r in pipe[d]]
        for i, (label, c, m) in enumerate(STAGES):
            ax.plot(x, [r[2 + i] for r in pipe[d]], color=c, marker=m, label=label)
        vals = [v for r in pipe[d] for v in r[2:6]]
        log_axes(ax, x, min(vals), max(vals))
        node_labels(ax, pipe[d])
        ax.set_title(f"({tag}) {d} by stage", loc="left")
        ax.set_xlabel("MPI ranks"); ax.set_ylabel("time (s)")

    # (f) alignment imbalance: slowest rank over the average rank
    for d in DATASETS:
        if not rank[d]:
            continue
        c, m = DATASET_STYLE[d]
        f.plot([r[1] for r in rank[d]], [r[5] for r in rank[d]], color=c, marker=m)
    xs = sorted({r[1] for d in DATASETS for r in rank[d]})
    f.set_xscale("log", base=2)
    f.xaxis.set_major_locator(FixedLocator(xs)); f.xaxis.set_minor_locator(NullLocator())
    f.xaxis.set_major_formatter(FuncFormatter(plain))
    f.set_xlim(xs[0] / 1.15, xs[-1] * 1.15); f.set_ylim(0.7, 7)
    f.axhline(1, color=GREY, lw=1.0, ls=":")
    f.set_title("(f) Alignment imbalance", loc="left")
    f.set_xlabel("MPI ranks"); f.set_ylabel("slowest / average rank")
    # dataset legend for (a) and (f); it sits in (f), where there is room
    f.legend(*a.get_legend_handles_labels(), loc="upper left", title="datasets in (a), (f)",
             title_fontsize=7, handlelength=1.6, labelspacing=0.3)

    h, l = panels[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.0), columnspacing=1.5)
    fig.tight_layout(rect=(0, 0, 1, 0.955), w_pad=1.2, h_pad=1.2)
    fig.savefig(os.path.join(HERE, "strong_scaling.pdf"))
    plt.close(fig)


def vs_fastga():
    pipe = by_dataset(read_table("table_pipeline.tex"))
    fga = by_dataset(read_table("table_fastga.tex"))      # threads, 3 steps, total
    fig, axes = plt.subplots(1, 4, figsize=(7.0, 2.25))
    for ax, d in zip(axes, DATASETS):
        ft, fy = [r[0] for r in fga[d]], [r[4] for r in fga[d]]
        px, py = [r[1] for r in pipe[d]], [r[6] for r in pipe[d]]
        ax.plot(ft, fy, color=GREY, ls="--", marker="s", label="FastGA (threads)")
        ax.plot(px, py, color=BLUE, marker="o", label="DistributedPanGA (ranks)")
        log_axes(ax, sorted(set(ft + px)), min(fy + py), max(fy + py))
        ax.set_title(d, loc="left")
        ax.set_xlabel("cores (threads or ranks)")
    axes[0].set_ylabel("total time (s)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.92), w_pad=0.8)
    fig.savefig(os.path.join(HERE, "vs_fastga.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    strong_scaling()
    vs_fastga()
