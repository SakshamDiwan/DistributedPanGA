#!/usr/bin/env python3
"""Figures and table numbers for the Results section.

Reads the per-run tables written by jobs/report/run_one.sh and fastga_one.sh:
  results.tsv  one row per pga-mpi-merge run
  fastga.tsv   one row per FastGA run
(default: ~/pga-report-results/, the copy that is not purged).

Writes into this directory:
  strong_scaling.pdf   total time vs ranks, and primates14 time per stage
  vs_fastga.pdf        pipeline vs FastGA at equal core counts

Each configuration is represented by its median run: the run whose total
time is the median of the repeats (lower middle for an even count). All
values shown for that configuration come from that one run, so its stage
times add up to its total. Figures show the same values as the tables.
Rerun after new jobs finish:  python3 report/fig/make_figs.py
"""
import csv
import os
import statistics
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/pga-report-results")

# Categorical slots 1-3 (fixed order) + a distinct marker per dataset.
DATASETS = [  # name in tsv, label, color, marker
    ("yeast7", "yeast7", "#2a78d6", "o"),
    ("primates4", "primates4.chr6", "#eb6834", "s"),
    ("primates14", "primates14.chr6", "#1baf7a", "^"),
]
STAGES = [  # column, label, color, marker
    ("extract_s", "extraction", "#2a78d6", "o"),
    ("sort_s", "distributed sort", "#eb6834", "s"),
    ("merge_s", "seed generation + redistribution", "#1baf7a", "^"),
    ("alignments_s", "alignment + PAF output", "#eda100", "D"),
]
INK, MUTED, GRID = "#222222", "#6b6b6b", "#e3e3e3"

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "text.color": INK, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "lines.linewidth": 1.5, "lines.markersize": 5,
    "legend.frameon": False, "pdf.fonttype": 42,
})


from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator, ScalarFormatter  # noqa: E402


def plain(ax, xs):
    """log axes with plain-number ticks: x at the measured points, y at 1-2-5."""
    ax.xaxis.set_major_locator(FixedLocator(sorted(set(xs))))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    lo, hi = ax.get_ylim()
    ys = []
    for mant in [(1, 2, 5), (1, 2, 3, 5), (1, 1.5, 2, 3, 4, 5, 7), (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8)]:
        ys = [m * 10 ** e for e in range(-2, 5) for m in mant if lo <= m * 10 ** e <= hi]
        if len(ys) >= 3:
            break
    ax.yaxis.set_major_locator(FixedLocator(ys))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def med(rows, col):
    """median, min, max of a column over rows (ignores NA)."""
    v = [float(r[col]) for r in rows if r[col] not in ("NA", "")]
    return statistics.median(v), min(v), max(v)


def median_run(rows):
    """the run whose total_s is the median (lower middle if even)."""
    ok = sorted((r for r in rows if r["total_s"] not in ("NA", "")),
                key=lambda r: float(r["total_s"]))
    return ok[(len(ok) - 1) // 2]


def run_dir(row):
    return os.path.join(SRC, row["dataset"], f"n{row['nodes']}-r{row['ranks']}",
                        f"rep{row['rep']}")


def fix_sort_wait(row):
    """Count the wait for the slowest extractor as extraction, not sort.

    Each rank stops its extraction timer when it finishes, then waits in the
    sort's first collective (MPI_Allreduce) for the slowest extractor. The
    per-stage maxima therefore count that wait twice: in extraction (slowest
    rank) and in sort (fastest rank). From the per-rank [o1-0] lines:
        sort = max_r(extract_r + sort_r) - max_r(extract_r)
    i.e. the time from the last rank finishing extraction to the end of the
    sort. Extraction keeps max_r(extract_r). Total is recomputed.
    """
    pairs = []
    with open(os.path.join(run_dir(row), "merge.err")) as f:
        for line in f:
            if line.startswith("[o1-0] rank"):
                kv = dict(x.split("=", 1) for x in line.split() if "=" in x)
                pairs.append((float(kv["extract_s"]), float(kv["sort_stage_s"])))
    if not pairs:
        return
    e_max = max(e for e, _ in pairs)
    row["sort_raw_s"] = row["sort_s"]
    row["sort_s"] = f"{max(e + s for e, s in pairs) - e_max:.6f}"
    row["total_s"] = f"{sum(float(row[c]) for c in ('extract_s', 'sort_s', 'merge_s', 'alignments_s')):.6f}"


runs = [r for r in read(os.path.join(SRC, "results.tsv")) if r["rc"] == "0"]
for r in runs:
    fix_sort_wait(r)
P = defaultdict(list)  # (dataset, nodes, ranks) -> rows
for r in runs:
    k = (r["dataset"], int(r["nodes"]), int(r["ranks"]))
    if k[0] == "primates14" and k[2] // k[1] != 16:
        continue            # report uses 16 ranks per node for primates14
    P[k].append(r)

fpath = os.path.join(SRC, "fastga.tsv")
F = defaultdict(list)  # (dataset, threads) -> rows
if os.path.exists(fpath):
    for r in read(fpath):
        if r["rc"] == "0":
            F[(r["dataset"], int(r["threads"]))].append(r)


def series(ds):
    """Main strong-scaling series: 1 node for yeast/primates4,
    16 ranks per node for primates14."""
    keys = sorted(k for k in P if k[0] == ds)
    if ds == "primates14":
        keys = [k for k in keys if k[2] // k[1] == 16]
    return keys


def errbar(ax, x, rows, col, color, marker, label=None, ls="-", fill=True):
    """plot the median run's value of col for each configuration (no error bars)."""
    m = [float(median_run(r)[col]) for r in rows]
    ax.plot(x, m, color=color, marker=marker, ls=ls, label=label,
            mfc=color if fill else "white", mec=color)
    return m


# ---------------------------------------------------------------- figure 1
fig, axs = plt.subplots(2, 2, figsize=(6.5, 5.0))
a = axs[0][0]
for ds, lab, col, mk in DATASETS:
    keys = series(ds)
    if not keys:
        continue
    x = [k[2] for k in keys]
    m = errbar(a, x, [P[k] for k in keys], "total_s", col, mk, lab)
    a.plot(x, [m[0] * x[0] / n for n in x], color=col, ls=":", lw=1)
a.set_xscale("log", base=2)
a.set_yscale("log")
a.set_xlabel("MPI ranks")
a.set_ylabel("total time (s)")
a.set_title("(a) Strong scaling (dotted: ideal)", loc="left")
plain(a, [1, 2, 4, 8, 16, 32, 64])
a.legend(loc="upper left")

# (b)-(d): time per stage for each dataset, one shared legend.
for b, letter, (ds, lab, _, _) in zip([axs[0][1], axs[1][0], axs[1][1]], "bcd", DATASETS):
    keys = series(ds)
    x = [k[2] for k in keys]
    for c, slab, col, mk in STAGES:
        errbar(b, x, [P[k] for k in keys], c, col, mk, slab)
    b.set_xscale("log", base=2)
    b.set_yscale("log")
    plain(b, x)
    if ds == "primates14":
        b.set_xticklabels([f"{n}\n({n // 16} node{'s' if n > 16 else ''})" for n in x])
    b.set_xlabel("MPI ranks")
    b.set_ylabel("time (s)")
    b.set_title(f"({letter}) {lab} by stage", loc="left")
h, l = axs[1][1].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=4, columnspacing=1.2, handlelength=1.6)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(os.path.join(HERE, "strong_scaling.pdf"))
if os.environ.get("PREVIEW"): fig.savefig(os.path.join(os.environ["PREVIEW"], "strong_scaling.png"), dpi=140)

# ---------------------------------------------------------------- figure 2
have = [d for d in DATASETS if any(k[0] == d[0] for k in F)]
if have:
    fig, axes = plt.subplots(1, len(have), figsize=(6.5, 2.4), squeeze=False)
    for ax, (ds, lab, _, _) in zip(axes[0], have):
        tk = sorted(k for k in F if k[0] == ds)
        errbar(ax, [k[1] for k in tk], [F[k] for k in tk], "total_s",
               "#6b6b6b", "s", "FastGA (threads)", ls="--")
        pk = series(ds)
        errbar(ax, [k[2] for k in pk], [P[k] for k in pk], "total_s",
               "#2a78d6", "o", "DistributedPanGA (ranks)")
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xlabel("cores (threads or ranks)")
        ax.set_ylabel("total time (s)")
        ax.set_title(lab, loc="left")
        plain(ax, [k[1] for k in tk] + [k[2] for k in pk])
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(os.path.join(HERE, "vs_fastga.pdf"))
    if os.environ.get("PREVIEW"): fig.savefig(os.path.join(os.environ["PREVIEW"], "vs_fastga.png"), dpi=140)

# ---------------------------------------------------------------- tables
# Table bodies only (rows); the table environment and caption live in the
# section file, so wording stays editable there.
LABEL = {d[0]: d[1] for d in DATASETS}
ORDER = {d[0]: i for i, d in enumerate(DATASETS)}


def fmt(v):
    return f"{v:,.2f}"


rows = []
prev = None
for k in sorted(P, key=lambda k: (ORDER.get(k[0], 9), k[2], k[1])):
    ds, nodes, ranks = k
    if prev is not None and ds != prev:
        rows.append("\\midrule\n")
    prev = ds
    mr = median_run(P[k])
    cells = [fmt(float(mr[c])) for c in
             ("extract_s", "sort_s", "merge_s", "alignments_s", "total_s")]
    rows.append(f"{LABEL.get(ds, ds)} & {nodes} & {ranks} & " + " & ".join(cells) + " \\\\\n")
with open(os.path.join(HERE, "table_pipeline.tex"), "w") as f:
    f.write(f"% Generated by make_figs.py from {SRC} -- do not edit.\n")
    f.writelines(rows)

rows = []
prev = None
for k in sorted(F, key=lambda k: (ORDER.get(k[0], 9), k[1])):
    ds, t = k
    if prev is not None and ds != prev:
        rows.append("\\midrule\n")
    prev = ds
    mr = median_run(F[k])
    cells = [fmt(float(mr[c])) for c in
             ("gixmake_s", "fastga_s", "alntopaf_s", "total_s")]
    rows.append(f"{LABEL.get(ds, ds)} & {t} & " + " & ".join(cells) + " \\\\\n")
with open(os.path.join(HERE, "table_fastga.tex"), "w") as f:
    f.write(f"% Generated by make_figs.py from {SRC} -- do not edit.\n")
    f.writelines(rows)

def per_rank_align(row):
    """alignment time (search_n_s + search_c_s) of every rank in one run,
    read from that run's merge.err ([pr1] lines)."""
    d = os.path.join(SRC, row["dataset"], f"n{row['nodes']}-r{row['ranks']}",
                     f"rep{row['rep']}", "merge.err")
    t = []
    with open(d) as f:
        for line in f:
            if line.startswith("[pr1] rank"):
                kv = dict(x.split("=", 1) for x in line.split() if "=" in x)
                t.append(float(kv["search_n_s"]) + float(kv["search_c_s"]))
    return t


rows = []
prev = None
for k in sorted(P, key=lambda k: (ORDER.get(k[0], 9), k[2], k[1])):
    ds, nodes, ranks = k
    if ds == "yeast7":
        continue
    if prev is not None and ds != prev:
        rows.append("\\midrule\n")
    prev = ds
    t = per_rank_align(median_run(P[k]))
    mean = sum(t) / len(t)
    rows.append(f"{LABEL.get(ds, ds)} & {nodes} & {ranks} & {fmt(max(t))} & "
                f"{fmt(mean)} & {fmt(min(t))} & {max(t) / mean:.2f} \\\\\n")
with open(os.path.join(HERE, "table_perrank.tex"), "w") as f:
    f.write(f"% Generated by make_figs.py from {SRC} -- do not edit.\n")
    f.writelines(rows)

# Weak scaling: runs with ~equal input bases per rank, and the alignment
# output per rank. Aligned bases and alignment counts are dataset totals from
# FastGA's gold PAF (every run's output is identical to it), divided by ranks.
FIXTURES = os.path.join(os.path.dirname(os.path.dirname(HERE)), "tests", "fixtures")
GOLD = {"yeast7": "yeast/gold.paf", "primates4": "primates4/gold.paf",
        "primates14": "primates14/gold.paf"}
INPUT_BP = {"yeast7": 83_409_131, "primates4": 672_453_404, "primates14": 2_630_931_800}
WEAK = [("yeast7", 1, 1), ("primates4", 1, 8), ("primates14", 2, 32),
        ("yeast7", 1, 2), ("primates4", 1, 16), ("primates14", 4, 64)]
gold_tot = {}
for ds, rel in GOLD.items():
    n = bp = 0
    with open(os.path.join(FIXTURES, rel)) as f:
        for line in f:
            c = line.split("\t", 11)
            n += 1
            bp += int(c[10])          # PAF column 11: alignment block length
    gold_tot[ds] = (n, bp)
rows = []
for i, k in enumerate(WEAK):
    if k not in P:
        continue
    if i == 3:
        rows.append("\\midrule\n")
    ds, nodes, ranks = k
    t = per_rank_align(median_run(P[k]))
    n, bp = gold_tot[ds]
    rows.append(f"{LABEL[ds]} & {ranks} & {INPUT_BP[ds] / ranks / 1e6:,.1f} & "
                f"{bp / ranks / 1e6:,.1f} & {n / ranks:,.0f} & "
                f"{fmt(sum(t) / len(t))} & {fmt(max(t))} \\\\\n")
with open(os.path.join(HERE, "table_weak.tex"), "w") as f:
    f.write(f"% Generated by make_figs.py from {SRC} and the gold PAFs -- do not edit.\n")
    f.writelines(rows)

print(f"wrote figures and tables from {SRC} ({len(runs)} pipeline runs, "
      f"{sum(len(v) for v in F.values())} FastGA runs)")
