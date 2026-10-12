# O1-0 results — per-rank extraction profile and pre-sort shuffle (yeast7)

- Job 59679320 (2026-10-11, nid004266, 1 node, exclusive).
- Commit `abaea3f` (branch `o1/0-instrumentation`); tree clean apart from the
  untracked fixture symlinks.
- Same configuration as the baseline (`baseline.md`): `-c 16 --cpu-bind=cores`,
  `OMP_NUM_THREADS=1`, `cc`, `-DLCPs`.
- Raw data: `$SCRATCH/pga-o1-0-59679320/` (per run: `merge.err` `[o1-0]` lines,
  `time.rank*`). Log: `~/pga-jobs/o1/logs/pga-o1-0.59679320.log`.

Each finding below is labelled **measured** or **hypothesis**. A hypothesis is
consistent with the data but has not been isolated.

## Correctness of the instrumented build

- `test-mpi` at ws 1/2/3/4/8 passes.
- `verify-phase3-stage2` on yeast7 is ALL PASS (Tier 1 + Tier 2) at N=1/2/4/8.
- 18 of 18 e2e runs `gold=MATCH`: 43,126 alignments and 0 duplicate lines in
  every run.
- Login node: `make test` (42,471 checks, 0 failures), `make test-asan`, and
  small-fixture `verify-phase1/2` all pass.
- Wall times agree with job 59676625 within 1–2%, so the instrumentation adds
  no measurable overhead.

## 1. Work distribution: balanced (measured)

| N | records per rank (min–max) | contigs per rank | extract_s per rank (min–max, 3 reps) |
|---|---|---|---|
| 1 | 64.38M | 119 | 49.8–50.1 |
| 2 | 31.84–32.54M | 58–61 | 30.1–31.1 |
| 4 | 15.83–16.31M | 29–31 | 23.2–24.3 |
| 8 | 7.01–8.82M | 13–16 | 28.2–33.5 |

- The whole-contig split is within ±11% at N=8 (`compute_work_split`). Extract
  time per rank tracks records per rank: rank 7 has the fewest records
  (7.01M) and is the fastest (29.0 s).
- **Imbalance does not explain the N=8 regression.** Even the lightest rank
  takes longer at N=8 than every rank takes at N=4.

## 2. Cost per record grows with N at fixed work per record (measured)

| N | mean extract_s | µs per record (aggregate) | aggregate throughput (M records/s on the node) |
|---|---|---|---|
| 1 | 49.95 | 0.78 | 1.29 |
| 2 | 30.58 | 0.95 | 2.11 |
| 4 | 23.73 | 1.48 | 2.71 |
| 8 | 31.97 | 3.97 | 2.01 |

Node-level extraction throughput tops out at about 2–2.7 M records/s,
whatever the rank count.

## 3. Kernel time does not shrink with per-rank work (measured)

From `/usr/bin/time -v` per rank. This covers the whole process, not just
extraction. Reproducible to ±0.5 s over 3 reps.

| N | user s / rank | system s / rank | system s, all ranks |
|---|---|---|---|
| 1 | 89.6 | 26.4 | 26 |
| 2 | 46.5 | 20.7 | 41 |
| 4 | 24.2 | 18.8 | 75 |
| 8 | 13.4 | 29.8 | 238 |

User time scales as ~1/N. System time per rank stays at 19–30 s while each
rank's work falls 8×, so total kernel time grows 9× for the same total work.

## 4. One `lseek` syscall per extracted record (measured, small fixture)

- `strace -f -c ./pga-extract small6` (38,023 records, serial extraction
  path), run with the `.bps` on both GPFS home and Lustre scratch:
  - **38,035 `lseek`**, 38,026 of them `SEEK_SET` on the `.bps` fd.
  - Only **84 `read`**.
- Each `pack_forward_kmer`/`pack_rc_kmer` call reaches `Get_Contig_Piece`,
  whose `fseeko(b, off, SEEK_SET)` (`lib/GDB.c:1908`) issues an `lseek`. The
  bytes themselves are then served from the stdio buffer (Lustre
  `st_blksize` = 4 MiB, so there are almost no `read` calls).
- `pack_reads == records` exactly on every rank at every N, so yeast7 at N=1
  makes ~64.4M `lseek`s.
- The per-call count is measured on `pga-extract` (the tuple path). That
  `extract_to_records` behaves the same is inferred from shared code
  (`scan_contig` and the same `pack_*` functions). It is not straced.

## 5. Stage 2 at N=8 (measured)

- The `sort_s` increase at N=8 (1.7 → 5.1 s in the baseline) is
  **Allreduce wait**: `hist_s` reaches 4.5–4.9 s on the ranks that finish
  extraction first.
- At N=8, `exchange_s` is 0.03–0.04 s, `local_sort_s` 0.22–0.25 s and
  `pack_s` 0.05–0.07 s. At N=1 they are 0.12 s, 1.61 s and 0.39 s.

## 6. Pre-sort shuffle (evidence, scoped)

- Setup: `PGA_PRESORT_SHUFFLE` ∈ {1,2,3} at N=1 and N=4 (6 runs).
- Every run is `gold=MATCH`, and the **raw, unsorted per-rank PAF files are
  byte-identical** to the unshuffled run's.
- Scope: yeast7, two rank counts, three seeds. This is evidence that output
  does not depend on the pre-sort order of equal k-mers on this input. It is
  not a proof for other inputs, which is why every later step keeps the PAF
  `cmp` gate.
- Cost: the shuffle itself takes ~3.5 s at N=1 (random access over 1 GB).

## Hypotheses (not yet isolated)

- **H1: the per-record `lseek` syscall is the dominant extraction cost at
  N=1.**
  - Support: 26.4 s system time at N=1 ≈ 64.4M syscalls × ~0.4 µs.
  - Open: some of that system time may come from later stages (merge/align
    `/tmp` I/O, `ALNtoPAF`).
  - Isolate by: comparing system time and `extract_s` with O1a, which removes
    these calls, or by `strace -c` on one rank's extraction only.
- **H2: the N=8 regression is contention among processes making `lseek`s on
  the same file**, e.g. a per-inode lock in the Lustre client.
  - Support: system time per rank is flat while user time scales; aggregate
    throughput saturates; work is balanced.
  - Not excluded:
    - node-level effects (CPU frequency, shared caches), though these would
      inflate user time, which scales cleanly;
    - contention in other kernel paths.
  - Isolate by: the same runs with the GDB copied to node-local `/tmp`
    (tmpfs), and with per-rank private copies (inode-level vs client-level).
    O1a removes the syscalls whatever the cause.
- **H3: after O1a, extraction becomes a small, compute-bound share of the
  time.** No evidence yet. Measure before scoping O1b.

## Consequences for prioritisation

1. **O1a first.** It directly removes the measured per-record `lseek`, and both
   H1 and H2 predict a large drop in extraction time. Its job should also
   include the H2 control (baseline binary with the GDB on `/tmp`) so the
   cause is attributed rather than assumed.
2. **O1b is re-scoped after O1a measurements.**
   - If extraction is still a large share, threads help.
   - If H2 holds, threads must not reintroduce per-record syscalls on the
     shared file; O1a's window reads are one per contig/window.
3. **O1c:** the shuffle evidence allows the simpler bucket-major pack order
   (with the PAF gate kept).
   - Motivation is mainly memory: peak RSS at N=1 is 4.4 GiB, and one record
     copy is ~1 GiB.
   - Runtime gain is limited: `hist_s` + `pack_s` at N=1 are ~0.45 s.
4. **O1d stays deferred.** The measured exchange is ≤ 0.12 s per rank on one
   node, so there is nothing to overlap at this scale.
