# Decision gate 1 — should we build HySortK-style static task dispatch?

Data:
- Step 3 job 59715282: PR 1 + PR 2, worktree commit `385cb08`.
- PR 1 job 59712760.
- Raw data: `$SCRATCH/pga-step3-l1-59715282/`; per run `balance.txt` and
  `pairs.rank*.tsv`.

## Correctness (all runs)

| dataset | gold | ranks | result |
|---|---|---|---|
| scerevisiae8 | 50,046 | 1, 2, 4, 8 | MATCH ×4 (N=2 aborted before PR 2) |
| primates4.chr6 | 1,547,485 | 4, 8, 16 | MATCH ×3 (N=4 needed chunked rounds: 2 stage-2 and 4 seed rounds) |

## Measured

**Wall time.** primates4.chr6 takes 179.4 / 98.8 / 71.6 s at N = 4 / 8 / 16.
The alignment stage is 80–83% of it (CSV `alignments_s` 149.1 / 80.0 / 60.1 s).

**Alignment-time imbalance** (`search_seeds`, max/mean over ranks):

| dataset | N=2 | N=4 | N=8 | N=16 |
|---|---|---|---|---|
| yeast7 | 1.05 | 1.07 | 1.11 | — |
| scerevisiae8 | 1.02 | 1.05 | 1.10 | — |
| primates4.chr6 | — | 1.14 | 1.26 | **1.89** |

Imbalance grows with N on every dataset. **This part of the hypothesis holds.**

**Seeds as the cost of a task** (one task = one `align_contigs` call):

| dataset | R² (linear) | Spearman | R² (log–log) |
|---|---|---|---|
| yeast7 | 0.89 | 0.67 | — |
| scerevisiae8 | 0.54 | 0.69 | — |
| primates4.chr6 | **0.19** | 0.93 | 0.63 (time ∝ seeds^1.15) |

- On primates4 the most expensive task has 75,487 seeds and takes 23.4 s.
  That is 5.2% of all task time. Tasks with ~40 M seeds take ~5 s.
- Half of all task time is in 147 of 14,290 tasks.
- **Seed count does not identify the expensive tasks. This part of the
  hypothesis fails** on the dataset where balancing matters.

## Model (computed from the measured task times; not a run)

Slowest rank's alignment time. Units are unordered contig pairs (today's
co-location unit). Each run's own task times are used.

| primates4 | today | static greedy by seeds (HySortK-style) | dynamic queue, seeds-descending order | perfect static (true times) | largest single task |
|---|---|---|---|---|---|
| N=4 | 123.5 s | 118.5 | 108.5 | 108.5 | 22.9 |
| N=8 | 69.1 s | **70.5 (worse)** | 62.9 | 54.5 | 22.8 |
| N=16 | 53.4 s | 48.3 | 43.9 | 28.1 | 23.4 |
| N=32 (extrapolated) | — | — | 33.2 | 14.0 → floor 23.4 | 23.4 |

## Verdict

- **Do not build static seed-count dispatch** (the planned PR 3). It recovers
  little (−4% at N=4, −10% at N=16) and is *worse* than today at N=8.
- HySortK's dispatch works because a sorting task's size *is* its cost. For
  alignment tasks it is not.
- The imbalance is real and growing, so the scheduling problem stays worth
  solving. It needs one of the following:
  1. **A better cost estimate.** Find out what makes the 23 s task slow, e.g.
     long high-identity or repeat alignments. A cheap feature known before
     alignment (seed diagonal concentration, contig identity) might predict
     it.
  2. **Dynamic scheduling** (a task queue or work stealing): about −9% at N=8
     and −18% at N=16 in the model. It still depends on ordering.
  3. **Splitting the largest tasks.** Beyond ~16 ranks the single largest
     task (~23 s) is the floor for any scheduler, so further scaling needs
     tasks smaller than one contig pair. That must keep alignments identical
     to FastGA, which is the hard part.
- **Next (diagnostic only, no change to the assignment):**
  - characterise the expensive tasks: contigs, seeds per diagonal, aligned
    bases from the PAF;
  - run primates14.chr6 to see whether the pathology generalises.
  - Then choose between 1 and 2, with 3 as the scaling extension.
