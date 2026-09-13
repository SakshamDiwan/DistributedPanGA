---
name: pga-validate
description: Run the DistributedPanGA phase verifiers and interpret their results. Use when asked to validate, verify, regression-test, or check correctness of extraction, the local or distributed sort, or LCP values, or when asked about gold fixtures and their provenance.
---

Full reference, including prerequisites, per-check semantics and gold
provenance: `docs/agent/validation.md`. Allocation details:
`docs/agent/perlmutter.md`.

## Workflow

1. Check prerequisites first — a missing fixture is the common failure. All
   five verifiers need `tests/fixtures/yeast/` artifacts, which `.gitignore`
   excludes and Git will never supply. The three MPI verifiers also need an
   active allocation.
2. Run, from the repository root, cheapest first:
   - `make verify-phase1` — extraction vs `gold.tuples`. Serial.
   - `make verify-phase2` — extraction + `msd_sort` vs `gold.tuples`. Serial.
   - `make verify-phase3-stage1` — gather-to-rank-0 baseline. Needs `srun`.
   - `make verify-phase3-stage2` — distributed sort. Needs `srun`.
   - `make verify-phase3-stage2-lcp` — LCP regression. Needs `srun`, and runs
     `make clean` first, so keep it separate from any other build or run.
3. Scope with `RANK_COUNTS` (default `1 2 4 8`), e.g.
   `RANK_COUNTS="1 2" make verify-phase3-stage2`. It has **no effect** on
   `scripts/runtime_test_phase4b.sh`.
4. For stage 2, require **both** checks to pass: Tier 2 (k-mer order across the
   unsorted rank-ordered concatenation) proves the sort was really distributed,
   and Tier 1 (multiset equality after sorting vs gold) proves nothing was lost
   or duplicated. Neither implies the other.
5. For LCP, remember `PGA_DUMP_LCP=1` switches `pga-mpi` to the 5-column dump
   (LCP `$4`, k-mer `$5`); the oracle is `gold.lcp.kmers.simple`, built by
   `scripts/build_lcp_oracle_simple.sh`, and GIXshow's `lcp` column is **not**
   a valid oracle.
6. End-to-end runs are the separate `bash scripts/runtime_test_phase4b.sh` —
   inside an allocation, never under `srun`.
7. Report which targets ran at which rank counts with the PASS/FAIL lines.
   Never describe a source-derived command as executed, and never cite the
   historical 43,126-alignment yeast7 figure as reproduced.
