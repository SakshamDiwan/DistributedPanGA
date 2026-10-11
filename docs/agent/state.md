# State / handoff

Last updated: 2026-10-11 (after O1-0).

## Current objective

Week 1 is complete: verified commit, yeast7 baseline and O1 design. Next is
O1, revised in `o1-design.md` (revision 2) so that the HySortK contribution is
the destination-aware staging (O1c), not just threading. The first code step
is **O1-0** (instrumentation and the tie-order experiment), then **O1a**.

## Commit and environment

- `b4f5d3b` (merged `main`), branch `chore/week1-baseline`.
- Build and run under `module load cpu; module unload darshan` with `cc`
  (`build.md`).
- FastGA tools are in `~/FASTGA` (`178afcb`).
- yeast7 fixtures live in `$SCRATCH/pga-fixtures/`, symlinked from
  `tests/fixtures/`. Scratch is purged after ~8 weeks idle; regenerate with
  `~/pga-jobs/week1/pga_w1_B.sbatch`.
- Job scripts and logs: `~/pga-jobs/week1/`. The user submits all Slurm jobs.

## Validation executed (all green)

- **Small fixture** (job 58873003): unit (42,471 checks), ASan/UBSan, MPI tier
  at ws 1–8, all phase verifiers, e2e.
- **yeast7:** phase1/2/stage1/stage2/LCP at N=1/2/4/8 (jobs 58873975 and
  59676624).
- **yeast7 e2e:** sorted, duplicate-preserving PAF `cmp` vs FastGA gold
  matches in all 16 runs (jobs 58873975 and 59676625).
- Details: `tests/RESULTS.md`. Timings: `baseline.md`.

## Key baseline facts

- Wall time: 113.3 / 67.1 / 42.9 / 44.4 s at N = 1/2/4/8.
- `extract_s`: 48.8 / 31.3 / 24.1 / **33.9** s. Extraction regresses at N=8
  and is 76% of wall time there; every other stage scales.

## Open unknowns

- Is the cause of the extraction regression per-seed `Get_Contig_Piece` I/O,
  Lustre contention, or work-split imbalance? O1-0 answers this.
- Does the PAF depend on the pre-sort order of equal k-mers (`msd_sort` is not
  stable)? O1-0 runs the shuffle experiment, which decides O1c's ordering.
- HySortK has no licence file: reuse ideas only.

## Follow-up defects (recorded, not fixed)

- **Fixed since 2026-09-10:** items 1 (send-displacement guard) and 2 (R-01
  seam LCP).
- **Still open:**
  - `int8_t` world-size cap of 127.
  - Two `ncontig²` maps.
  - 128 MiB adapter index per rank.
  - Non-fatal `ALNtoPAF` failure.
  - Dead `cleanup` label and `read_file_all`.
  - `O(ixlen)` `GoTo_Kmer_Index`.
  - Script portability.
  - Dead doc references.
- **New:**
  - `Read_GDB` returns 1 on error, but drivers test `< 0`, so the check can
    never fire.
  - The CSV reports only `MPI_MAX` per stage, so `total_s` overstates wall
    time at N=8.

## O1-0 (done, job 59679320)

Findings are in `o1-0-results.md`.
- **Measured:**
  - Work is balanced across ranks.
  - Each extracted record costs one `lseek` syscall (`fseeko` in
    `Get_Contig_Piece`).
  - Per-rank system time stays at 19–30 s whatever N is.
  - The N=8 `sort_s` rise is Allreduce wait.
  - The shuffled pre-sort order gives an identical PAF (yeast7, N=1/4,
    3 seeds).
- **Hypotheses:** the `lseek`s dominate extraction (H1), and contention on
  the shared file causes the N=8 regression (H2).

## Next action

Implement **O1a** (pack k-mers from the scan window) on a branch off
`o1/0-instrumentation`. Its measurement job must include the H2 control:
baseline and O1a binaries, with the GDB on Lustre and on node-local `/tmp`,
at N=1 and N=8. Gates as in `o1-design.md` §5.
