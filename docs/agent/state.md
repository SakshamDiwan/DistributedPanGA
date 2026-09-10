# State / handoff

Last updated: 2026-09-10.

## Current objective

Agent context for the DistributedPanGA checkout is set up. The next piece of
work is the optimization proposal (O1–O5 in `research-plan.md`), which cannot
start from measurements yet because no gold fixtures exist here.

## Inspected commit

`ce64406b45908451bd161466ebab0646cc84b92d` ("completed implementation with
phase wise validation"), on branch `chore/agent-context`. That branch was
created from this commit; `main` and `origin/main` point at the same commit, so
the checkout matched the grounding revision exactly and no reconciliation was
needed. The working tree was clean before this documentation change.

Tracked top-level content: `Makefile`, `.gitignore`, `src/`, `lib/`,
`scripts/`, `notes/`. No tracked README or tests upstream. Several script
comments reference `implementation_plan.md` and `src/phase4b_runtime_testing_plan.md`,
which **do not exist** in the repository (`*.md` is `.gitignore`d), so those
cross-references are dead ends.

## Validation actually performed

Executed on a Perlmutter login node:

| Check | Result |
|---|---|
| `make all` | exit 0; all six binaries produced |
| `make -B CC=cc pga-mpi` | exit 0; 409128 B vs 388888 B from `mpicc` (different MPI linked) |
| `make -n LCPS=0 pga-mpi`, `make -n CFLAGS=-O2 pga-mpi` | both "up to date" — staleness pitfall confirmed |
| `make -Bn CFLAGS=-O2 pga-mpi` | `-DLCPs`, `-I` and `-fno-strict-aliasing` all absent — override pitfall confirmed |
| driver-only compile of `src/pga-mpi-merge.c` | two dead-code warnings (below) |
| MPI-4 large-count probe | `MPI_Alltoallv_c`/`MPI_Alltoall_c` link; `MPI_Count` is 8 B; MPI 4.1 via `cc`, 4.0 via `mpicc` |
| Slurm/account probe | `m4341` has `interactive` and `debug` QOS; CPU node = 128 physical / 256 logical cores, ~503 GiB |

**Not executed:** every `make verify-phase*` target, `scripts/make_gold.sh`,
`scripts/runtime_test_phase4b.sh`, and any `srun`. No job was submitted and no
allocation was held. Everything in `validation.md` is source-derived unless it
says otherwise.

Binaries from the `make all` above are present in the repository root and are
`.gitignore`d. They were built with `mpicc` except `pga-mpi`, which was last
rebuilt with `cc` — **do not mix them in one measurement**; run
`make clean && make -B CC=cc all` before benchmarking.

## Missing local prerequisites

These block all correctness and performance work, in priority order:

1. `tests/fixtures/yeast/` does not exist — no `yeast7.1gdb`, `yeast7.gix`,
   `gold.tuples`, `gold.lcp.kmers.simple`, or `gold.paf`. `tests/` is
   `.gitignore`d, so Git will never provide them.
2. `yeast7.fa.gz`, the input `scripts/make_gold.sh` needs, is also absent.
3. No FastGA tool installation: `~/FASTGA` does not exist, and neither does
   `/global/u2/s/sdiwan/FASTGA` (hardcoded in `runtime_test_phase4b.sh`).
   `FAtoGDB`, `GIXmake`, `GIXshow`, `FastGA` and `ALNtoPAF` are all off `PATH`.
   `ALNtoPAF` is required at runtime by `pga-mpi-merge`.
4. `notes/GIXmake_dump` is a tracked prebuilt binary; rebuilding it needs
   `~/FASTGA` sources.

## Open unknowns

- Whether the historical yeast7 gold set is recoverable, or whether gold must
  be regenerated from `yeast7.fa.gz` with a freshly built FastGA. The reference
  figure of **43,126 alignments** is historical and has not been reproduced
  here; regenerating gold with a different FastGA revision may not reproduce
  it, which would be a provenance question rather than a bug.
- Which FastGA revision produced the existing gold artifacts. Nothing in the
  repository records it.
- Whether `cc` (cray-mpich 9.1.0) or `mpicc` (MPICH 9.0.1) was used for the
  historical timings. This changes the linked MPI library and so the numbers.

## Follow-up defects found during inspection

Recorded, not fixed. Roughly by severity.

1. **Unguarded cumulative send displacement** —
   `src/pga-mpi-merge.c:174` (`redistribute_seeds`). Per-destination send sizes
   and the cumulative *receive* total are checked against `INT_MAX`, but
   `send_displs[r] = (int) cum` casts a running `int64` with no cumulative
   *send* guard. Many destinations each under 2 GiB but summing above it will
   overflow silently instead of aborting. Add the symmetric guard.
2. **Wrong seam LCP past an empty rank** — `src/dsort.c:253`
   (`verify_cross_rank_boundaries`). A rank with zero records sends
   `my_last = 0x00...` to its successor. The successor's boundary comparison
   then passes vacuously, and more importantly the LCP seam fixup at
   `dsort.c:434` computes `simple_kmer_lcp_bases(0x00..., my_first)` against
   that all-zero sentinel rather than against the true global predecessor,
   which lives further back. The neighbour exchange only ever looks one rank
   away. Yields a wrong LCP byte on one record per empty-rank seam; needs a
   test with an empty intervening rank.
3. **`world_size` capped at 127** — `src/contig_assignment.c:82`, `int8_t`
   storage. Rejected cleanly with a message, so it is a limit rather than a
   bug, but it caps a single job at 127 ranks. Note the `N*N` allocation at
   line 75 happens *before* that check.
4. **Quadratic contig-map storage, twice over** —
   `contig_assignment.c:75` allocates `ncontig²` bytes for the upper triangle,
   and `pga-mpi-merge.c:406` allocates a second `NCONTS²` symmetric copy. Both
   are live simultaneously.
5. **128 MiB per-rank adapter index** — `src/kmer_adapter.c:39`, `ibyte = 3`
   gives `2^24` `int64` entries, allocated before the record table and merge
   cache. Fixed cost regardless of slice size, so it dominates at high rank
   counts.
6. **Converter failure is not fatal** — `src/pga-mpi-merge.c:519`. `ALNtoPAF`
   failure prints a message and `main` still returns 0, so a zero exit can
   accompany missing or truncated PAF output.
7. **Dead code in the driver** — `src/pga-mpi-merge.c:555` label `cleanup` is
   defined but never jumped to (so the error paths that presumably should reach
   it currently `MPI_Abort` instead, leaking nothing but bypassing cleanup),
   and `read_file_all` at line 102 is unused. Both warn under `-Wall`.
8. **`O(ixlen)` `GoTo_Kmer_Index`** — `src/kmer_adapter.c:143`, a 2^24 linear
   scan. Harmless at `NTHREADS = 1`; becomes a real cost if merge threading is
   ever raised, as the comment there notes.
9. **Script portability** — `make_gold.sh` assumes `$HOME/DistributedPanGA` and
   `~/FASTGA`; `runtime_test_phase4b.sh` hardcodes another user's FastGA path,
   truncates the timing CSV and deletes `out.rank*.paf` on every run, and uses
   fixed `/tmp` filenames. See `perlmutter.md`.
10. **Dead documentation references** — scripts cite `implementation_plan.md`
    and `src/phase4b_runtime_testing_plan.md`, neither of which exists.

## Next action

Smallest useful next step, in order:

1. Obtain or rebuild the FastGA tools and put them on `PATH`; confirm with
   `command -v ALNtoPAF`.
2. Restore or regenerate `tests/fixtures/yeast/` and record the provenance
   fields listed in `validation.md`.
3. In one interactive allocation
   (`salloc --account=m4341 --constraint=cpu --qos=interactive --nodes=1 --time=00:30:00`),
   run `make verify-phase1` and `make verify-phase2` — serial, no `srun`, and
   the cheapest possible confirmation that extraction and the local sort still
   match gold. Then `verify-phase3-stage2` at `RANK_COUNTS="1 2"`, and
   `verify-phase3-stage2-lcp` last because it runs `make clean` first.
4. Only then take a baseline measurement, recording every field in
   `research-plan.md`.

Do not begin O1–O5 implementation before step 3 passes: without gold data
there is no oracle to detect a regression.
