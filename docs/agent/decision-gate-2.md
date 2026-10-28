# Decision gate 2 — what predicts alignment cost, and what scheduling can achieve

- Jobs 59717508 (1 node) and 59717520 (2 nodes, the first multi-node run).
- Binary `addf993` (`diag/align-cost`: PR 1 + PR 2 + `PGA_SCRATCH_DIR` +
  per-task features).
- Raw data: `$SCRATCH/pga-diag-p14-n16-59717508/` and
  `$SCRATCH/pga-p14-n32-59717520/`.

## Correctness (measured)

| dataset | ranks / nodes | wall | gold | result |
|---|---|---|---|---|
| primates4.chr6 | 8 / 1 | 100.9 s | 1,547,485 | MATCH |
| primates4.chr6 | 16 / 1 | 75.0 s | 1,547,485 | MATCH |
| primates14.chr6 | 16 / 1 | 1,440.9 s | 4,577,556 | MATCH (4 stage-2 and 6 seed rounds) |
| primates14.chr6 | 32 / 2 | 1,221.9 s | 4,577,556 | MATCH (4 and 4 rounds) |

- **Peak RSS per rank** (`overall_peak`): primates14 11.5 GB at N=16 and
  8.7 GB at N=32.
- **The alignment stage is 96–97% of primates14 wall time** (CSV
  `alignments_s` 1,389 s at N=16 and 1,188 s at N=32).
- **Doubling ranks from 16 to 32 gives only 1.18×.**

## Imbalance (measured, `search_seeds` max/mean)

primates4: 1.27 (N=8), 1.89 (N=16). primates14: **1.95 (N=16), 3.35 (N=32)**.

## What predicts a task's time (per-task fits; same tasks at every N)

| feature | known | primates4 R² / R²(log) / Spearman | primates14 (N=32) R² / R²(log) / Spearman |
|---|---|---|---|
| seeds | before | 0.19 / 0.62 / 0.93 | 0.06 / 0.62 / 0.94 |
| distinct diagonals | before | 0.21 / 0.55 / 0.90 | 0.16 / 0.54 / 0.90 |
| max seeds on one diagonal | before | 0.07 / 0.36 / 0.90 | 0.04 / 0.48 / 0.90 |
| contig length product | before | 0.25 / 0.44 / 0.84 | 0.04 / 0.40 / 0.78 |
| log-space model on all of the above | before | — / 0.63 / 0.94 | — / 0.63 / — |
| bases covered by non-redundant alignments | after | 0.45 / **0.93** / 0.93 | 0.71 / **0.91** / 0.95 |

Cost follows **output size** (aligned bases), which is not known in advance.
Nothing measurable beforehand explains more than R²(log) ≈ 0.63.

## Where the cost is

- primates4:
  - The three most expensive tasks involve two small contigs (885 kbp and
    379 kbp).
  - One pair's non-redundant alignments cover 152 Mbp. That is
    tandem-repeat behaviour: few seeds, because FREQ=10 drops frequent
    k-mers, but enormous alignment work.
- primates14:
  - **Self-alignment tasks (contig vs itself) are 49% of all alignment work.**
  - The largest single task is contig 0 against itself: 773–788 s, 7% of
    all work. Its alignments cover 7.4 Gbp of a 14 Mbp contig.
  - Total work divided by the largest task is ≈ 14, so the alignment stage
    cannot use much more than ~14 ranks without splitting tasks.

## Scheduling models (from measured task times; not runs)

| slowest rank (s) | measured today | static by seeds | dynamic by seeds | **dynamic, self-first by contig length** | perfect | largest task |
|---|---|---|---|---|---|---|
| primates4 N=8 | 69 | 76.3 | 55.2 | **54.7** | 54.7 | 22.8 |
| primates4 N=16 | 53 | 50.8 | 40.4 | **31.0** | 28.4 | 23.7 |
| primates14 N=16 | 1,351 | 1,312 | 959 | **876** | 788 | 788 |
| primates14 N=32 | 1,159 | 1,246 | 826 | **774** | 773 | 773 |

- "Dynamic" means each rank takes the next task from one shared order as it
  becomes free.
- The self-first order uses only facts known before alignment: whether the
  task is a self-alignment, and the contig length.

## Conclusions

1. **Static size-based dispatch (the HySortK mechanism) is ruled out**
   (confirms gate 1). It is no better, and sometimes worse, than today's
   length-based assignment, because size does not predict cost.
2. **Dynamic scheduling with a self-first order is close to optimal.** It
   cuts the slowest rank's alignment time by 33–42% at 16–32 ranks in the
   model. It is the next thing to build. Prediction to verify: the measured
   wall time drops correspondingly.
3. **Beyond ~16 ranks, the largest self-alignment task is the floor for any
   scheduler on primates14.** Scaling further needs **splitting
   self-alignment tasks** while keeping alignments identical to FastGA. That
   is the open research problem: it needs a study of how `align_contigs`
   handles a contig against itself (diagonal bands, redundancy elimination)
   before any design.
4. **Relation to HySortK:** what transfers is the task abstraction (many
   independent tasks, more than ranks). The static, size-based dispatcher
   does not. The evidence for that is itself a reportable result.
