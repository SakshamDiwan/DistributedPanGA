# Research plan — proposal mapping and test priorities

Source commit inspected: `ce64406`. This file maps the optimization proposal
onto concrete code locations and records what must be measured. It does **not**
authorize implementing any of it; the current task is documentation only.

## Baseline to preserve

- `FREQ = 10` is the baseline and stays. The "roughly 10× haplotype count"
  heuristic is a research choice for larger datasets, not an instruction to
  change the baseline.
- `KMER=40`, `TMER=12`, `SMER=8`, `NUM_BUCK=1024` are fixed; changing any of
  them invalidates every gold artifact.
- `NTHREADS = 1` in the merge/align drivers and `nthreads = 4` in
  `sort_records` are the current settings. Report them with every measurement.
- Alignment and chaining behaviour (`entwine`, `align_contigs`, `search_seeds`,
  `pair_sort_search`, `la_merge` in `src/fastga_pipeline.c`) is treated as
  correct-as-is. Changes to surrounding I/O must preserve it, and "the process
  exited 0" is not evidence that it did.

## O1 — HySortK extraction integration

Target: `extract_to_records` / `scan_contig` (`src/extract.c:65`), currently
serial within each rank and the extract stage of the timing CSV.

Constraints any replacement must satisfy:

- Preserve **syncmer selection** as implemented, including the `mz == min4`
  fall-through tie case and the left-end rescan branch. A different minimizer
  scheme changes the seed set and therefore the alignments.
- Preserve the **occurrence payload**: every emitted seed carries
  `(contig, position, strand)`. HySortK-style pipelines that produce aggregate
  per-k-mer *counts* are not a substitute — the merge stage needs individual
  occurrences to form pairs, and `recompute_buck` needs per-contig counts.
- Preserve **LCP**: the downstream adapter copies each record's LCP byte into
  the table (`src/kmer_adapter.c:68`) and `new_self_merge_thread` reads it
  (`suf1[CBYTE]`, `p[-2]`) to bound adaptamer groups. Dropping or changing LCP
  semantics silently changes seeding.
- Preserve both emissions per position `j` (forward stores `j`, RC stores
  `j+TMER`) and the exact coordinate convention.

Acceptance gate: `make verify-phase1` and `make verify-phase2` byte-identical
against `gold.tuples`, then `verify-phase3-stage2` and
`verify-phase3-stage2-lcp`.

## O2 — large counts, within-contig splitting, bucket balance

Three separable pieces.

**(a) Count/displacement widths.** Two independent places:

- `src/dsort.c:organize_send_buffer` guards both the per-rank byte count and
  the running total against `INT_MAX` and aborts — correct but a hard ceiling.
  `run_stage2` guards the cumulative receive total likewise.
- `src/pga-mpi-merge.c:redistribute_seeds` guards each **per-destination** send
  size (`send_sizes64[r] > INT_MAX`) and the cumulative **receive** total, but
  computes `send_displs[r] = (int) cum` from a running `int64` **with no
  equivalent guard on the cumulative send total**. With many destinations each
  under 2 GiB but summing over it, displacements overflow silently rather than
  aborting. Fix the guard before anything else here.
- The `world_size == 1` shortcut in `redistribute_seeds` skips MPI entirely, so
  it avoids the limit on one rank only. It does not remove any multi-rank
  limit, and it means the 1-rank path is not exercising the exchange code at
  all — a correctness blind spot as well as a scaling one.

MPI-4 large-count availability, **verified locally on 2026-09-10**: a test
translation unit using `MPI_Alltoallv_c` / `MPI_Alltoall_c` compiles and links
under both wrappers, `sizeof(MPI_Count) == 8`, and `MPI_VERSION`/`SUBVERSION`
report **4.1 via `cc`** (cray-mpich 9.1.0) and **4.0 via `mpicc`** (MPICH
9.0.1). So the `_c` large-count variants are a real option here, not
speculative — but they are only worth adopting together with `int64`
count/displacement arrays throughout, and the alternative (chunking) keeps
portability. Re-verify on any environment change before committing to it.

**(b) Within-contig splitting.** `compute_work_split` (`src/work_split.c`)
assigns whole contigs, so extraction load is bounded below by the largest
contig and a single-contig genome cannot be split at all. Splitting inside a
contig must handle the scan's `KMER`/`TMER` boundary windows — a naive split
loses or duplicates seeds within `KMER` bases of each cut.

**(c) Bucket balance.** `record_bucket` uses only the first 5 bases → 1024
buckets, and `compute_destination_map` assigns contiguous bucket ranges by
global count. With `world_size` approaching 1024, or with skewed 5-base
prefixes, granularity itself becomes the limit. Any change must stay
deterministic and communication-free, since correctness depends on all ranks
computing the identical map.

## O3 — seed-density-aware assignment

`contig_assignment_build` currently weights pairs by `clen[i]*clen[j]`, a proxy
for work that ignores actual seed density.

The hard requirement is that **every rank must derive the identical map**.
Today that holds because the only inputs are replicated GDB metadata
(`gdb->ncontig`, `contigs[].clen`), so no communication is needed. Measured
seed densities are *not* replicated — each rank only sees its own slice after
redistribution. So a density-aware policy must first make the weights global,
e.g.:

- an `MPI_Allreduce` over a per-contig seed-count vector (length `ncontig`,
  `int64`) after the sort stage, then run the existing deterministic heap on
  the reduced vector; or
- assign on rank 0 and `MPI_Bcast` the map.

Either adds a synchronization point that does not exist today; state which and
account for its cost. Determinism also requires care with floating-point
weights — keep integers.

Also fix, or at least account for, the storage: `int8_t` caps `world_size` at
127, and two `ncontig²` byte maps are live at once
(`contig_assignment.c` upper-triangle map plus the symmetric `RankForPair` in
the driver). At 32k contigs that is ~1 GiB each.

## O4 — remove seed / `.las` / `.1aln` round trips

Current flow per rank: merge writes seeds to files under `/tmp`, they are read
back into memory for `MPI_Alltoallv`, written to a combined file, re-read by
`reimport_thread`, chained into `.las` scratch, merged into `.1aln`, and then
converted by an **external `ALNtoPAF` via `system()`**.

Opportunities and hazards:

- The file round trips are inside the `merge` timing bucket, so the CSV already
  attributes them; measure before assuming they dominate.
- Eliminating the external converter also removes the failure mode where
  `ALNtoPAF` is missing from `PATH` or fails and `main` still returns 0.
- `recompute_buck` exists only because redistribution invalidates the counts
  the merge wrote; an in-memory exchange could maintain them directly.
- Scratch is node-local `/tmp` with `TMPDIR` ignored, so any redesign must not
  assume cross-node file visibility.

## O5 — GPU

Deferred. Current resources are CPU-only and the account for this work is
**m4341** (CPU); `m4341_g` is out of scope. Do not add GPU code paths or CUDA
dependencies to the build.

## Benchmark record — required fields

For every reported measurement, record all of:

- Input dataset identity and a **hash** of the input and of the reference
  (`sha256sum` of `.fa.gz`, `.1gdb`, `gold.paf`).
- Source revision (`git rev-parse HEAD`) **and dirty state**
  (`git status --porcelain`); a dirty tree makes a number unciteable unless the
  diff is recorded.
- Compiler and MPI versions and which wrapper (`cc` vs `mpicc` — they link
  different MPI libraries here) and the exact expanded compile command
  (`make -Bn <target>`), including whether `-DLCPs` was present.
- Nodes, ranks, ranks-per-node, `--cpus-per-task`, `--cpu-bind`, and the
  effective per-rank thread counts (`nthreads=4` in the sort, `NTHREADS=1` in
  the merge).
- Slurm job ID.
- **Actual wall time measured separately** from the CSV's `total_s`, which is
  the sum of per-stage `MPI_MAX` values and overestimates end-to-end time.
- Memory: state the method (`sacct MaxRSS`, `/usr/bin/time -v`, or
  in-process instrumentation) — they are not interchangeable.
- Repetitions and dispersion, not a single run.
- Correctness outcome for the same configuration, and how it was checked.

Agreement across world sizes on a larger dataset is **consistency** evidence
only. It is not an independent oracle: a systematic error present at every
world size agrees with itself. An independent oracle means FastGA-derived gold
data.

## Future unit tests (do not implement now)

Priority order, all with **independently computed** expected values rather than
values captured from the current implementation:

1. **Packing and coordinates** — for a short synthetic contig, assert
   `pack_forward_kmer` bytes for `[j, j+40)` and `pack_rc_kmer` for
   `[j-28, j+12)` against hand-computed 2-bit packings; assert the stored
   positions are `j` and `j+TMER`; assert `-1` at both range boundaries.
2. **Semantic LCP** — `simple_kmer_lcp_bases` on hand-built pairs: identical
   (expect 40), differing in the first base (0), differing only in the last
   base (39), and a difference at each of the four positions within a byte.
3. **Record widths** — `compute_record_sizing` against hand-computed
   `post_bytes`/`cont_bytes` for boundary `maxctg`/`ncontig` values (exactly
   256, 65536, `2*ncontig` crossing a byte); assert the strand bit round-trips
   through the high bit of the last contig byte and does not corrupt the contig
   id.
4. **Whole-contig splits** — `compute_work_split` for `num_workers` greater than
   `ncontig`, one giant contig plus many small ones, and equal-length contigs;
   assert boundaries are non-decreasing, cover `[0, ncontig)` exactly, and that
   `post[num_workers] == seqtot`.
5. **Bucket maps** — `record_bucket` extracts exactly the first 5 bases;
   `compute_destination_map` is a total function onto `[0, world_size)`,
   produces non-decreasing `ksplit`, and behaves for empty buckets and
   `world_size > ` the number of non-empty buckets.
6. **Deterministic assignment** — `contig_assignment_build` for two
   equal-length contigs on two ranks must give `(0,0)→0, (0,1)→1, (1,1)→0`
   under the current tie policy (weight tie → `(i,j)` ascending → least-load
   heap with lower-rank tie-break). Verified by hand against the code, not by
   running it. Also assert identical output across repeated calls and that
   `world_size = 128` is rejected.

Distributed-path cases that deserve explicit tests once a harness exists:
empty and single-record rank slices, and an **empty intervening rank** (see the
seam defect in `state.md`).

Note that `tests/` is `.gitignore`d, so adding a test suite requires a
deliberate decision about how it is tracked.
