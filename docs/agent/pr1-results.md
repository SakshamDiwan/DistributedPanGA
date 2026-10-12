# PR 1 results — merge/alignment instrumentation and early frees

- Job 59712760 (2026-10-12, nid006677, 1 node, exclusive), commit `2e884ca`.
  The overall-peak/opt-out follow-up is `da9830d`.
- Baseline: `main` `fd4b6e9`. Configuration as `baseline.md`: `-c 16
  --cpu-bind=cores`, `NTHREADS=1`, `cc`.
- Raw data: `$SCRATCH/pga-pr1-59712760/`. Per run there is `merge.err` with
  `[pr1]`/`[pr1-mem]` lines, `pairs.rank*.tsv` and `balance.txt` (from
  `scripts/analyze_align_balance.py`).

## Correctness

- **Small fixture:** e2e at N=1/2/4/8 passes. The ASan/UBSan driver
  (`make pga-mpi-merge-asan`) at N=1/2/4 shows 0 ASan errors, 0 UBSan reports
  outside the vendored `lib/`, and gold MATCH.
- **yeast7** (gold 43,126): base N=1/8 and PR 1 N=1/2/4/8 ×3 all MATCH.
- **scerevisiae8** (gold 50,046): base N=1/8 and PR 1 N=1/4/8 ×3 all MATCH.
- **scerevisiae8 at N=2 aborts** on `main` and PR 1 alike: "outgoing N bytes
  exceed the int range of MPI_Alltoallv". That is the existing `INT_MAX`
  limit (~3.7 GB of seeds per rank), fixed by PR 2 (R-06).
- **Invariant:** for every run and rank, the per-task seeds sum to
  `seeds_n + seeds_c` (`analyze_align_balance.py` exits 1 otherwise).

## Memory (measured; peak RSS per rank, max over ranks)

| | `main` (`/usr/bin/time`) | PR 1 (`overall_peak`) | change | stage holding the peak |
|---|---|---|---|---|
| yeast7 N=1 | 4,463 MiB | 2,961 MiB | −34% | stage 2 (exchange) |
| yeast7 N=8 | 1,071 MiB | 721 MiB | −33% | seed redistribution |
| scerevisiae8 N=1 | 6,113 MiB | 3,791 MiB | −38% | alignment (sort array) |
| scerevisiae8 N=8 | 1,597 MiB | 1,195 MiB | −25% | seed redistribution |

- **Prediction vs measurement:** I predicted −45% at N=1 (~2.3 GiB). The
  measured gain is smaller because once the alignment-stage buffers are gone,
  stage 2 (local + send + receive buffers) becomes the peak.
- **At N>1** the peak is the old `redistribute_seeds`, which holds full send
  and receive buffers. PR 2 streams that path.
- **Method caveat:**
  - PR 1's per-stage reset of VmHWM also resets `ru_maxrss`, so
    `/usr/bin/time` read 18–41 MiB for PR 1 runs. That number is wrong.
  - Use `overall_peak_kb` instead, or set `PGA_MEM_STAGE_RESET=0`.
- Wall time is unchanged within noise: yeast7 N=1 65.68 → 65.83 s, N=8
  10.56 → 10.57 s.

## Where the alignment-stage time goes (measured, N=8, rep 1)

| per rank | yeast7 mean / max (s) | scerevisiae8 mean / max (s) |
|---|---|---|
| pair_sort_search total | 5.82 / 6.63 | 7.84 / 8.70 |
| — align (`search_seeds`, N+C) | 2.96 / 3.29 | 3.96 / 4.49 |
| — seed sort (`rmsd_sort`) | 1.85 / 2.19 | 2.51 / 3.19 |
| — seed reimport from `/tmp` | 1.00 / 1.18 | 1.35 / 1.72 |
| merge core (`self_adaptamer_merge`) | 1.27 / 1.29 | 1.67 / 1.70 |
| seed redistribution + `recompute_buck` | 0.85 / 1.03 | 1.24 / 1.54 |
| `ALNtoPAF` | 0.05 | 0.05 |

Rows that combine N and C phases give the sum of the per-phase means and of
the per-phase maxima. Those maxima can come from different ranks, so the "max"
of a combined row is an upper bound.

The seed path (redistribute, recompute, reimport, sort) costs as much as the
alignment itself (~3.7 s of yeast7's 10.6 s wall time).

## Imbalance and predictability (measured)

Max/mean across ranks, rep 1:

| | seeds | align (`search`) | pair search | merge+align |
|---|---|---|---|---|
| yeast7 N=2 | 1.074 | 1.049 | 1.064 | 1.050 |
| yeast7 N=4 | 1.098 | 1.074 | 1.087 | 1.069 |
| yeast7 N=8 | 1.181 | 1.111 | 1.139 | 1.110 |
| scerevisiae8 N=4 | 1.089 | 1.058 | 1.075 | 1.058 |
| scerevisiae8 N=8 | 1.206 | 1.102 | 1.109 | 1.089 |

- **Per task** (one `align_contigs` call: type, anchor, partner), task
  seconds vs seeds:
  - yeast7: 28,155 tasks, R² = 0.89–0.90, Spearman 0.67.
  - scerevisiae8: 37,903 tasks, **R² = 0.54**, Spearman 0.69.
- Seeds are a fair predictor on yeast7 and a mediocre one on scerevisiae8,
  which contains the near-identical S288C and SGDref.
- **Per rank** (N=8), search time vs seeds R² = 0.66 (yeast7) and 0.51
  (scerevisiae8).

## Scheduling bound (model, not measurement)

- **Method:** longest-first assignment of the *measured* tasks to the same
  number of ranks.
- **Oracle** (true task times): max = 1.000× mean at every N. There are
  thousands of small tasks, so perfect balance is possible.
- **Seed-count estimate:** max = 1.002–1.036× mean.
- **At N=8 the seed-based bound would cut the slowest rank's align time:**
  - yeast7: 3.15 → 2.90 s;
  - scerevisiae8: 4.16 → 3.89 s.
- Balancing the whole pair search (sort and reimport also scale with seeds)
  bounds the gain at the gap between max and mean: about 0.8 s, **~6–8% of
  wall time at N=8**, with the imbalance growing with N.

These yeast-scale numbers alone do not clear the decision gate (they show
imbalance growing with N, but a modest absolute gain). The primate datasets
(Step 3) decide it.
