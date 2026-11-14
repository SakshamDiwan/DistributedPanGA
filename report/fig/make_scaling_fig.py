"""Cross-dataset scaling summary (report Section 5.4).

Values are copied from fig/table_pipeline.tex and fig/table_perrank.tex
(Saksham, make_figs.py) and from the human runs (Table tab:human, jobs
61255434 / 61244457; slowest/avg from the [pr1] search_n_s + search_c_s
per rank). Rerun: `module load python; python3 make_scaling_fig.py`.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

total = {  # ranks -> total s
    "yeast7": {1: 79.97, 2: 42.81, 4: 22.64, 8: 12.69},
    "primates4.chr6": {4: 211.47, 8: 114.02, 16: 80.76},
    "primates14.chr6": {16: 1463.94, 32: 1218.54, 64: 1229.09},
    "hsapiens90.chr6": {32: 291.17, 64: 180.07},
}
imbalance = {  # ranks -> slowest / average alignment time per rank
    "primates4.chr6": {4: 1.09, 8: 1.23, 16: 1.70},
    "primates14.chr6": {16: 1.93, 32: 3.23, 64: 6.59},
    "hsapiens90.chr6": {32: 1.32, 64: 1.49},
}
style = {
    "yeast7": ("#4d9078", "o"), "primates4.chr6": ("#f78154", "s"),
    "primates14.chr6": ("#b4436c", "^"), "hsapiens90.chr6": ("#253237", "D"),
}

fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.6))
for name, d in total.items():
    r = sorted(d)
    c, m = style[name]
    a.plot(r, [d[r[0]] / d[x] for x in r], color=c, marker=m, ms=4, label=name)
    a.plot(r, [x / r[0] for x in r], color=c, ls=":", lw=0.8)
a.set_xscale("log", base=2); a.set_yscale("log", base=2)
a.set_xlabel("Ranks", fontsize=8); a.set_ylabel("Speedup over smallest run", fontsize=8)
a.set_title("(a) Total time", fontsize=8)
for name, d in imbalance.items():
    r = sorted(d)
    c, m = style[name]
    b.plot(r, [d[x] for x in r], color=c, marker=m, ms=4, label=name)
b.set_xscale("log", base=2)
b.axhline(1, color="grey", lw=0.6, ls=":")
b.set_xlabel("Ranks", fontsize=8); b.set_ylabel("Slowest / average rank", fontsize=8)
b.set_title("(b) Alignment imbalance", fontsize=8)
for ax in (a, b):
    ax.tick_params(labelsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xticks([1, 2, 4, 8, 16, 32, 64]); ax.set_xticklabels(["1", "2", "4", "8", "16", "32", "64"])
b.set_xlim(3, 80)
a.set_yticks([1, 2, 4, 8]); a.set_yticklabels(["1", "2", "4", "8"])
h, l = a.get_legend_handles_labels()
fig.legend(h, l, fontsize=7, ncol=4, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.02))
fig.tight_layout()
fig.savefig("scaling_summary.pdf", bbox_inches="tight")
fig.savefig("/tmp/claude-112970/-global-u2-h-hsc53/d0a40387-3000-4ba0-809a-53bbb355352f/scratchpad/scaling_summary.png",
            bbox_inches="tight", dpi=110)
