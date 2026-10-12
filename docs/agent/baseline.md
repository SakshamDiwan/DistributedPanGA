# Week-1 baseline — yeast7, `pga-mpi-merge`, 1 CPU node

The reference point for every O1 measurement. Rerun with
`sbatch ~/pga-jobs/week1/pga_w1_D.sbatch` (unchanged) for comparison.

## Configuration

| Field | Value |
|---|---|
| Slurm job | 59676625 (2026-10-11, node nid004745, `regular`, `--exclusive`) |
| Source | `b4f5d3b`. Dirty only in docs, `tests/RESULTS.md` and `scripts/verify_phase3_stage2.sh`; no `src/` or `lib/` change |
| Binary | `pga-mpi-merge` sha256 `8f4a40cc…31bbaf` |
| Build | `module load cpu; module unload darshan; make -B CC=cc all` → `cc -O3 -Wall -Wno-unused-result -fno-strict-aliasing -I./lib -I./src -DLCPs` |
| Toolchain | gcc-native 14, PrgEnv-gnu 8.7.0, cray-mpich 9.1.0, cpe 26.03 |
| Node | 2 × AMD EPYC 7763 (128 cores / 256 hardware threads) |
| Placement | ranks N ∈ {1,2,4,8} on one node; `srun -n N -c 16 --cpu-bind=cores` |
| Threads | `OMP_NUM_THREADS=1`, merge `NTHREADS=1`, sort `nthreads=4` (compiled in) |
| Parameters | `FREQ=10`, `KMER=40`, `TMER=12`, `SMER=8`, `NUM_BUCK=1024`, `SCAN_MAX=10^7` |
| Input | `yeast7.fa.gz` sha256 `385292…7db20571c992b0ccfefb1a827db97e08` (7 strains × 16 chr, 83.5 Mbp, 119 GDB contigs) |
| Reference | `gold.paf` sha256 `cd4230…56be` (FastGA `178afcb`, 43,126 alignments). See `tests/fixtures/yeast/MANIFEST.txt` |
| Fixture storage | `$SCRATCH` (Lustre), through symlinks in `tests/fixtures/` |
| Repetitions | 3, interleaved (rep → N=1,2,4,8) |

## Results (mean ± sample sd, n = 3)

| N | wall (s) | extract_s | sort_s | merge_s | align_s | CSV total_s | max rank RSS (MiB) |
|---|---|---|---|---|---|---|---|
| 1 | 113.26 ± 0.30 | 48.83 ± 0.14 | 2.55 ± 0.01 | 14.20 ± 0.01 | 46.74 ± 0.03 | 112.32 ± 0.15 | 4460 |
| 2 | 67.11 ± 0.57 | 31.31 ± 0.64 | 2.29 ± 0.13 | 8.57 ± 0.01 | 24.86 ± 0.05 | 67.02 ± 0.74 | 3647 |
| 4 | 42.87 ± 0.15 | 24.10 ± 0.14 | 1.67 ± 0.24 | 4.40 ± 0.00 | 12.75 ± 0.06 | 42.93 ± 0.33 | 1910 |
| 8 | 44.43 ± 0.13 | 33.93 ± 0.05 | 5.06 ± 0.05 | 2.33 ± 0.01 | 6.65 ± 0.02 | 47.96 ± 0.13 | 1072 |

**Correctness:** all 12 runs exit 0 with 43,126 alignments and 0 duplicate
lines, and every rank writes a non-empty file with no `ALNtoPAF` failures. The
sorted, duplicate-preserving aggregate PAF is byte-identical to FastGA gold
(`gold=MATCH`) in every run.

### How each column is measured

- **wall:** `date` around `srun` in `pga_w1_e2e_lib.sh`, so it includes launch
  and teardown. It is the end-to-end number to cite.
- **Stage columns:** `pga-mpi-merge.timing.csv`. Each value is the `MPI_MAX`
  over ranks (`pga-mpi-merge.c:524`), so per-rank spread is invisible.
  - `extract_s` covers `record_buffer_init` + `extract_to_records`.
  - `sort_s` covers all of `run_stage2` (histogram, exchange, local sort, LCP,
    boundary checks).
  - `total_s` is the sum of the stage maxima, which is why at N=8 it exceeds
    wall time.
- **RSS:** the per-rank `/usr/bin/time -v` "Maximum resident set size", maximum
  over ranks. `sacct MaxRSS` is not used. It disagrees badly: 7.3 GB at N=1,
  where it includes the wrapper and samples, and ~120 MB at N=8, where the
  sampler misses the peak.

## Observations (inputs to O1)

- **Extraction stops scaling and then regresses:** 48.8 → 31.3 → 24.1 → 33.9 s.
  At N=8 it is 76% of wall time and slower than N=2, even though each rank
  scans ~1/8 of the bases.
  - The other stages scale near-linearly: merge 14.2 → 2.3 s, align
    46.7 → 6.65 s.
  - Candidate causes:
    - per-seed `Get_Contig_Piece` seek+read in `pack_*_kmer`, with N ranks
      contending on one Lustre file;
    - whole-contig work-split imbalance.
  - O1-0 instrumentation must separate these: per-rank extract time, records
    per rank and read-call counts.
- **`sort_s` rises at N=8** (1.7 → 5.1 s). This is the first sign of
  exchange/boundary cost, or imbalance absorbed inside `run_stage2`. O1-0
  splits it into `exchange_s` and `local_sort_s`.
- **N=1 extraction costs ≈0.76 µs per record** (≈64.4M records).
- Run-to-run variance is small (sd ≤ 1.3% of the mean for wall), so 3 reps are
  enough to resolve changes of a few percent.
