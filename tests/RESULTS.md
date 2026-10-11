# Test results

**Branch:** `main` (merged PRs #3–#6)  **Commit:** `b4f5d3b`  **Tree:** clean source
(untracked docs only)
**Run:** 2026-09-25, Perlmutter, account `m4341`. Unit suite on a login node;
MPI tier, verifiers and e2e in batch job **58873003** (1 CPU node, `debug`).
Build: `module load cpu; module unload darshan; make clean && make -B CC=cc all`
(gcc 14.3.0, cray-mpich 9.1.0, `-DLCPs`). See `docs/agent/build.md` for why
darshan must be unloaded.

R-01, R-03 and R-04 are fixed on this commit and **no expected-failure markers
remain**.

## Unit suite — `make CC=cc test`

**42,471 checks, 0 failures, 0 expected failures**, exit 0, ~23 s.

| Binary | Covers | Checks |
|---|---|---|
| `test_record_split` | E-12, E-13, R-03 sizing | 112 |
| `test_extract` | E-01…E-11, R-03 extraction | 7,519 |
| `test_extract_window` | E-11 (`SCAN_MAX=64`) | 19 |
| `test_lcp_sort` | S-01…S-05b, R-04 | 4,028 |
| `test_dsort` | S-06…S-09 | 26,668 |
| `test_assign` | A-01…A-05 | 76 |
| `test_adapter` | K-01…K-05 | 4,034 |
| `test_exchange_plan` | R-02 | 15 |

`make CC=cc test-asan` (ASan + UBSan) clean: `test_lcp_sort --slow` 4,032
checks (includes S-05b, which forces a production buffer growth),
`test_extract` 7,519.

## MPI tier — `make test-mpi`, world sizes 1 2 3 4 8

| Mode | Result |
|---|---|
| `uniform` | pass at 1, 2, 3, 4, 8 (11 checks each) |
| `skew` | pass at 1, 2, 3, 4, 8 (the former ws=2 R-01 failure is gone) |
| `allempty` | pass at 1, 2, 3, 4, 8 |
| `seam` (ws=4) | **pass**, R-01: observed LCP 4, correct 4; slices `r0=4 r1=0 r2=4 r3=0` |
| M-06 consistency | ws 2, 3, 4, 8 all match ws=1 |

## Phase verifiers — small synthetic fixture

Against independent FastGA gold (see `MANIFEST.txt` for provenance):

| Check | Result |
|---|---|
| `verify-phase1` | PASS — extraction byte-identical to `gold.tuples` (38,023 tuples) |
| `verify-phase2` | PASS — extraction + `msd_sort` matches after lex sort |
| `verify-phase3-stage1` | ALL PASS at N ∈ {1, 2, 4, 8} |
| `verify-phase3-stage2` | ALL PASS at N ∈ {1, 2, 4, 8}, Tier 1 + Tier 2 |
| `verify-phase3-stage2-lcp` | ALL PASS at N ∈ {1, 2, 4, 8}; 36,128 distinct-k-mer LCPs match the oracle |

## End-to-end — `scripts/verify_e2e_small.sh`

`pga-mpi-merge` at world sizes 1, 2, 4 and 8; aggregate PAF compared against the
FastGA reference, sorted but not deduplicated.

| World size | Rank files | Aggregate | Result |
|---|---|---|---|
| 1 | 1 | 24 alignments | matches FastGA gold |
| 2 | 2 | 24 alignments | matches FastGA gold, identical to ws=1 |
| 4 | 4 | 24 alignments | matches FastGA gold, identical to ws=1 |
| 8 | 8 | 24 alignments | matches FastGA gold, identical to ws=1 |

The reference containing 24 alignments is a property of the fixture, not
evidence about this pipeline; the evidence is that the aggregate **matches** it
at every rank count.

## yeast7 — phase verifiers and end-to-end vs independent FastGA gold

Fixture and gold: `tests/fixtures/yeast/MANIFEST.txt`. Gold was regenerated in
job 58873948 with FastGA `178afcb` and reproduces the historical **43,126**
alignments. A contig-order variant gives a byte-identical sorted PAF. Same
commit, same `cc` binaries as above.

| Check | Job | Result |
|---|---|---|
| `verify-phase1` | 58873975 | PASS (64,378,260 tuples) |
| `verify-phase2` | 58873975 | PASS |
| `verify-phase3-stage1` | 58873975 | ALL PASS at N ∈ {1, 2, 4, 8} |
| `verify-phase3-stage2` | 59676624 | ALL PASS at N ∈ {1, 2, 4, 8}, Tier 1 + Tier 2 |
| `verify-phase3-stage2-lcp` | 58873975 | ALL PASS at N ∈ {1, 2, 4, 8}; 14,547,417 distinct-k-mer LCPs match the oracle |
| e2e: sorted, duplicate-preserving aggregate PAF `cmp` vs gold | 58873975, 59676625 | MATCH at N ∈ {1, 2, 4, 8}: 43,126 alignments, 0 duplicate lines; 16 of 16 runs |

`verify-phase3-stage2` first reported Tier 2 FAIL at line 342 for every N
(job 58873975). That was a verifier bug: awk compared the hex k-mers
numerically, so `…83e9` was read as 8.3e10. The script now forces a string
comparison under `LC_ALL=C`, and the rerun (59676624) passes. No pipeline code
changed.

Baseline timings: `docs/agent/baseline.md`.

## Defect status

| ID | Status |
|---|---|
| R-01 | fixed (PR #3); `seam` and `skew` ws=2 pass |
| R-02 | fixed; markers removed |
| R-03 | fixed (PR #4); markers removed |
| R-04 | fixed (PR #3); marker removed |

Small-fixture provenance re-checked: every sha256 in
`tests/fixtures/generated/MANIFEST.txt` matches, including the five FastGA
binaries.

## Notes

- MPI binaries built with `cc` under the default modules link darshan and
  exit 255 (`OneTextSchema … errno 95`); first attempt, job 58872495. Build and
  run with `module load cpu; module unload darshan`.

- Tests build with `$(CFLAGS)`, so `-DLCPs` matches the `LCPS=1` default.
  Without it `recalc_all_lcps` silently takes its 0/1 branch.
- MPI binaries need `MPICH_GPU_SUPPORT_ENABLED=0` in the default Perlmutter
  environment, which loads `craype-accel-nvidia80`; otherwise cray-mpich aborts
  with *"GPU_SUPPORT_ENABLED is requested, but GTL library is not linked"*.
  `tests/run_mpi.sh` and `scripts/verify_e2e_small.sh` set it.
- The fixture and gold are regenerated by `scripts/make_gold_small.sh`, verified
  reproducible: a from-scratch rebuild produced a byte-identical `gold.paf`.

See `tests/CATALOG.md` for what each test asserts.
