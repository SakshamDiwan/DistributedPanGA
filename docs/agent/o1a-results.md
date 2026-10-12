# O1a results — pack k-mers from the scan window (yeast7)

- Commit `f528988` (branch `o1/a-window-pack`) compared with the pre-O1a
  binary `abaea3f`. Both were built with `cc`, darshan unloaded, `-DLCPs`.
  - `abaea3f` binary: sha256 `fe998688…`, identical to the binary measured in
    O1-0.
  - O1a binary: sha256 `4390da26…`.
- Configuration is the same as `baseline.md`: 1 node, `-c 16 --cpu-bind=cores`,
  `OMP_NUM_THREADS=1`, GDB on Lustre (`$SCRATCH`).

## Correctness (all green)

| Gate | Where | Result |
|---|---|---|
| `make test` (42,476 checks, incl. new `o1a_window_pack_paths`), `make test-asan` (now including `test_extract_window` at `SCAN_MAX=64`) | login node | pass |
| small6 `verify-phase1/2`; `pga-extract` output byte-identical to pre-O1a | login node | pass |
| `test-mpi` ws 1/2/3/4/8 | job 59686567 | pass |
| yeast7 `verify-phase1` (byte-identical to `gold.tuples`), `verify-phase2` | job 59686567 | PASS |
| yeast7 `verify-phase3-stage1`, `verify-phase3-stage2` (Tier 1 + 2) at N=1/2/4/8 | job 59686567 | ALL PASS |
| yeast7 e2e sorted, duplicate-preserving PAF `cmp` vs gold | jobs 59686567 (2 runs) and 59691033 (18 runs) | MATCH in 20 of 20 |
| yeast7 `verify-phase3-stage2-lcp` | job 59691033 | ALL PASS at N=1/2/4/8; 14,547,417 distinct-k-mer LCPs match the oracle. Post-check rebuild is byte-identical to the measured binary |

## Performance (job 59691033, nid006079; mean ± sample sd, n = 3, interleaved)

| bin | N | wall (s) | extract_s | sort_s | merge_s | align_s | sys s (all ranks) | user s (all ranks) | max RSS/rank (MiB) |
|---|---|---|---|---|---|---|---|---|---|
| base | 1 | 113.80 ± 0.58 | 49.76 ± 0.46 | 2.54 | 14.15 | 46.52 | 26.7 | 88.9 | 4457 |
| o1a | 1 | **65.43 ± 0.05** | **1.48 ± 0.00** | 2.55 | 14.15 | 46.53 | 2.7 | 64.6 | 4464 |
| o1a | 2 | 36.29 ± 0.06 | 0.75 ± 0.01 | 1.45 | 8.56 | 24.80 | 5.3 | 68.4 | 3650 |
| o1a | 4 | 19.03 ± 0.06 | 0.38 ± 0.00 | 0.79 | 4.40 | 12.67 | 5.5 | 70.2 | 1906 |
| base | 8 | 39.11 ± 0.19 | 28.61 ± 0.09 | 4.47 | 2.45 | 6.68 | 199.3 | 108.6 | 1072 |
| o1a | 8 | **10.57 ± 0.03** | **0.21 ± 0.00** | 0.45 | 2.33 | 6.61 | 6.3 | 73.0 | 1070 |

### Measured

- **Extraction:** 33.6× faster at N=1 (49.76 → 1.48 s) and 136× faster at N=8
  (28.61 → 0.21 s).
- **Wall time:** 1.74× faster at N=1 and 3.70× faster at N=8.
- **Scaling:** O1a's strong scaling from 1 to 8 ranks is 6.2× in wall time
  (base: 2.9×). Extraction scales ~7×.
- **Window packing:** every record is packed from the window
  (`window_packs == records`, `pack_reads == 0` on every rank). All yeast7
  contigs fit in one `SCAN_MAX` window, so the fallback never fires.
- **Kernel time:** 26.7 → 2.7 s at N=1 and 199.3 → 6.3 s at N=8 (all ranks).
- **Allreduce wait gone:** the N=8 `sort_s` rise (4.47 → 0.45 s) disappears.
  It was wait time behind the slowest extractor.
- **Other stages unchanged:** merge, align and peak RSS are identical within
  noise. O1a does not touch those stages.
- **Node-to-node variation:** base at N=8 measured 39.1 s wall here against
  44.4 s in the baseline job (different node). Comparisons are therefore made
  only within one job, with interleaved repetitions.

### Hypotheses revisited

- **H1 (the per-record `lseek` dominated extraction): strongly supported.**
  Removing exactly those calls, and changing nothing else in the scan, removes
  97% of the N=1 extraction time.
- **H2 (N=8 regression from contention on the shared file): the regression is
  gone with the syscalls.** Kernel time per rank at N=8 is now ~0.8 s.
  - The precise kernel mechanism (Lustre client vs other) remains
    **unattributed**. The `/tmp` control did not complete (see below).
  - One uncontrolled data point: base, N=1, GDB on `/tmp` (job 59686567, a
    different node) took `extract_s` = 36.7 s against ~50 s on Lustre.
  - The question is moot for performance now, so it is not pursued further.

## Incident: `ALNtoPAF` hang with the GDB on `/tmp`

- **What happened (job 59686567):** base binary, N=1, with the GDB copied to
  node-local `/tmp`. The run completed every pipeline stage, including the
  alignment search (`43126 non-redundant aln's`), then hung inside the
  external `ALNtoPAF` call.
  - `out.rank0.paf` was left at 0 bytes.
  - It ran > 2 h until the job's time limit, and node nid006480 was left in
    `CG`.
- **Unknown:** the driver discards `ALNtoPAF` stderr (`2>/dev/null`) and gives
  it no timeout, so the cause is unknown.
- **Follow-up** (extends defect 6 in `state.md`): give the converter a timeout,
  keep its stderr, and make its failure fatal.
- The `/tmp` control is now opt-in (`RUN_TMP_CONTROL=1`) in
  `~/pga-jobs/o1/pga_o1_a2.sbatch`.

## Consequences for the roadmap

- **Extraction share of wall time:** 2% at N=8 (0.21 of 10.57 s) and 2% at
  N=1. The remaining time is `align_s` (pair_sort_search + `ALNtoPAF`, 63%)
  and `merge_s` (22%). Both are outside O1, and the alignment algorithms are
  frozen.
- **O1b (threaded extraction): deprioritised.** It can save at most 0.21 s at
  N=8 or 1.5 s at N=1 on yeast7. Revisit only for larger inputs, where
  per-rank extraction grows.
- **O1c (destination-aware staging):** the runtime upside on yeast7 is small
  (`hist_s` + `pack_s` ≈ 0.45 s at N=1). Its case is **memory** (peak 4.4 GiB
  per rank at N=1; each record copy ≈ 1 GiB) and scaling to larger
  pangenomes. Pursue it only with a memory-focused measurement.
- **O1d:** still deferred. The measured exchange is ≤ 0.13 s.
