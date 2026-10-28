# O1 design — HySortK-informed extraction and send-path redesign

Status: **complete (O1 closed 2026-10-12).** See `o1-final-report.md` for what was built, what was dropped and why. Revision 2 (2026-10-11), written
against `b4f5d3b`. It replaces revision 1 (2026-09-25), which treated O1 as
"thread the scan" and called that HySortK integration. Implementation must not
start until the week-1 gate below is green.

## 0. Gate and evidence

**Gate: green as of 2026-10-11.**

- Job 59676624: `verify-phase3-stage2` ALL PASS (Tier 1 + Tier 2) on yeast7
  after the Tier-2 awk fix. The earlier "failure" was awk comparing hex k-mers
  as numbers.
- Job 59676625 (baseline, N ∈ {1,2,4,8} × 3 reps): all 12 runs `gold=MATCH`.
  Table in `docs/agent/baseline.md`, which supersedes the single-run table
  below.

Already green: unit, ASan, MPI tier, small-fixture verifiers and e2e (job
58873003). Yeast phase1/2/stage1/LCP, plus the sorted, duplicate-preserving
PAF `cmp` against gold at N = 1/2/4/8 (job 58873975).

**Why extraction.** Job 58873975 timing CSV (one run each, so indicative only):

| N | extract_s | sort_s | merge_s | align_s | total_s |
|---|---|---|---|---|---|
| 1 | 48.8 | 2.5 | 14.2 | 46.5 | 112.0 |
| 2 | 36.1 | 2.2 | 8.5 | 24.8 | 71.7 |
| 4 | 20.2 | 1.3 | 4.4 | 12.6 | 38.4 |
| 8 | 29.0 | 4.3 | 2.3 | 6.6 | 42.3 |

- Extraction is the only stage that stops scaling: N=8 is slower than N=4, and
  extraction is ~70% of the total at N=8.
- At N=1 it costs ~0.76 µs per record (≈64.4M records). That is far above the
  scan arithmetic, which points at the per-seed `Get_Contig_Piece`
  seek+read in `pack_*_kmer`. This is a **hypothesis**; O1-0 measures it.
- `sort_s`, which includes the whole exchange, is 1–4 s. Communication is *not*
  the bottleneck at this scale, and that limits what O1d can buy (§5).

## 1. What HySortK actually does (commit `1542fe3`, `icpp2024`)

References are to `src/kmerops.cpp` unless noted.

- **Pipeline** (`kmer_count`, `src/hysortk.cpp`): `prepare_supermer` →
  `exchange_supermer` → `filter_kmer`. Parsing does *not* overlap sending.
- **Destination finding** (`FindKmerDestinationsParallel`, ~1010):
  - OpenMP `schedule(static)` over reads, one contiguous block of reads per
    thread.
  - `ParallelData::destinations[tid][read][kmer]` stores one int task-id per
    k-mer.
  - Destination = `murmur(canonical minimizer) % tot_tasks`. There is no
    counting pass, only `reserve`.
- **Tasks ≠ ranks:** `tot_tasks = avg_tasks · nprocs` (≈ 3 tasks per worker,
  oversubscribed).
  - Tasks are assigned to ranks *after* a global size reduction, by
    `BalancedDispatcher` (~1214) or `RoundRobinDispatcher`.
  - `HeavyHitterClassifier` turns oversized tasks into locally pre-counted
    lists. It is disabled when `EXTENSION=1`.
- **Per-task × per-thread output** (`ScatteredSupermers`, ~253):
  `supermers_[tid]` and `lengths_[tid]` per task. Each thread appends only to
  its own slot, so no locks are needed. Send order is thread 0..T-1, which is
  global read order because of the static split.
- **Exchange** (`TaskManager`, `SendingGroup`/`ReceivingGroup`, ~587–1007):
  - One blocking `MPI_Alltoallv` of per-task item counts.
  - Then rounds of `MPI_Ialltoall` with fixed 80 KB per-destination slots,
    double-buffered (4 × 80 KB × P per rank). The next batch is packed while
    one collective is in flight. Progress relies on `MPI_Wait` only, with no
    `MPI_Test` loop.
  - Each slot's last byte is a done flag, ANDed globally, so all ranks leave
    the loop together.
  - It is two-stage: stage 1 sends lengths, stage 2 sends bases. The receiver
    therefore preallocates exactly and needs no framing.
  - The `int` count is always 80000, so there are no `INT_MAX` limits.
- **Supermers** (`SupermerEncoder::encode`, ~1109):
  - A supermer is a maximal run of consecutive k-mers with the same *task id*,
    2-bit packed, capped at `MAX_SUPERMER_LEN` = 250 bases.
  - With `EXTENSION=1` it carries one `(pos, rid)` per supermer, and positions
    are rebuilt as `pos+i`.
  - The receiver re-expands directly into a preallocated `kmer_` array inside
    the receive loop (`receive_from_buffer_stage2`, ~484).
- **Local sort:** PARADIS (in place) or RADULS (out of place), chosen by free
  memory (`sort_decision`).

**Incompatible with our semantics (not adopted):**

- **Minimizer/murmur routing, `BalancedDispatcher`, heavy hitters.**
  - Our ownership must be *contiguous k-mer prefix ranges*. That is what makes
    the concatenation of rank outputs globally sorted (Tier 2), and it is what
    lets the merge see every occurrence of a k-mer on one rank. Routing by
    minimizer hash breaks both.
  - HySortK's task→rank map also depends on thread count, which breaks
    rank-count-invariant determinism.
- **Canonical k-mers.** They drop the strand; we emit forward and RC
  occurrences separately with distinct coordinates.
- **Counting/filtering.** The merge needs every occurrence, not counts.
- **PARADIS/RADULS.** Sorting is out of scope (no algorithm changes).

## 2. Generic technique vs HySortK architecture

| Step | Technique | Classification |
|---|---|---|
| O1a | Pack k-mers from the resident scan buffer instead of re-reading the file | Generic I/O hygiene. The HySortK parallel is narrow: `copy_bits` encodes from the already-resident 2-bit `DnaSeq` and never rereads input. |
| O1b | pthread pool over contigs, per-thread GDB reader | Generic multithreading. Its determinism scheme (§4) is ours, not HySortK's; HySortK gets determinism from a static split. |
| O1c | Destination-aware staging: work unit × routing-bucket organisation, routing decided after a global reduction, fused histogram, direct pack into per-destination send slots | **HySortK's core architecture**: `ParallelData`/`ScatteredSupermers` + count reduction + dispatch, applied with *our* deterministic bucket map in place of their dispatcher. |
| O1d | Fixed-slot batched `MPI_Ialltoall`, double buffering, done flags, two-stage exchange, receiver-side expansion | **HySortK's exchange architecture** (`SendingGroup`/`ReceivingGroup`). |
| — | Supermer compression | **Not applicable** (§6); a weak analogue is noted there. |

**Our 1024 prefix buckets play the role of HySortK's oversubscribed tasks.**
Both are fine-grained routing units whose mapping to ranks is decided after a
global reduction:

- Theirs: `MPI_Reduce` sizes → `BalancedDispatcher`.
- Ours: `MPI_Allreduce` histogram → `compute_destination_map`, which is
  deterministic and gives contiguous ranges because `select[]` is
  non-decreasing.

That correspondence is what lets O1c adopt their data organisation without
their routing.

## 3. Correctness invariants (all steps)

1. **Syncmer selection unchanged.** `scan_contig` (`src/extract.c:69-209`)
   keeps its control flow byte-for-byte: `mz==min4` falls through to emit
   (:176), left-end rescan (:164-173), `SCAN_MAX` windows with rolling state
   (:196-207).
2. **Every occurrence** `(kmer, contig, pos, strand)` is emitted through the
   existing record layout (`record_emit_cb`, :254-283; strand is the
   contig-field high bit).
3. **Coordinates.** Forward covers `[j, j+40)` and stores `j`, emitted iff
   `j ≤ len-40`. RC covers `[j-28, j+12)` and stores `j+12`, emitted iff
   `j ≥ 28`.
4. **Semantic LCP.** Extraction writes LCP byte 0. LCPs are produced only by
   `recalc_all_lcps` (`src/sort.c:38`) plus the seam fixup
   (`src/dsort.c:466-475`). Untouched.
5. **Pre-sort order.**
   - `msd_sort` is deterministic but **not stable**, so the arrangement of
     equal k-mers after sorting depends on the arrangement of its input.
   - Whether merge/PAF output depends on that is unverified.
   - Requirement until O1-0 answers it: **the byte stream each destination
     receives must be identical to today's** (source rank order, then
     extraction order within a source).
6. **Final alignments identical:** the sorted, duplicate-preserving aggregate
   PAF `cmp` against FastGA gold at N = 1/2/4/8.
7. **No changes** to sorting (`sort_records`, `msd_sort`), the adaptive merge,
   chaining or alignment.

## 4. Engineering concerns and how each is handled

- **Thread-safe GDB access.**
  - In EXTERNAL mode, `Get_Contig_Piece` does `ftello/fseeko/fread` on the
    shared `gdb->seqs` (`lib/GDB.c:1907-1914`), which races.
  - Fix: a per-thread reader. Shallow-copy the `GDB` struct with its own
    `FILE*` from `fopen(gdb->seqpath,"r")`. `Get_Contig_Piece` reads only
    `seqs`, `seqstate`, `ncontig` and `contigs`. Alternatively use
    `pread(fileno(...))`.
  - Never call `Load_Sequences`.
  - Test hook: the synthetic GDB has `seqpath == NULL` and an `fmemopen`
    stream (`tests/support/synth_gdb.c:109-113`). The reader factory therefore
    takes an "open another stream" callback that reopens `fmemopen(sg->bps)`
    in tests.
- **Other shared state.**
  - `Comp[]` (`src/tables.c:3`) is rewritten by every `run_scan` (:293).
    Initialise it once, before threads start.
  - `RecordBuffer` growth is per-thread or per-task, never shared.
  - `EPRINTF`/`Error_Buffer` and `EXIT` are not thread-safe. Error paths in
    threads set a flag; the main thread reports it.
  - `msd_sort` file statics: the sort stays one call per rank.
- **Determinism under dynamic scheduling.**
  - Work units are claimed dynamically (atomic counter, largest first), so
    *which thread* scans a unit is nondeterministic.
  - Therefore all output is keyed by **work unit**, never by thread. Each unit
    writes to its own segments, and assembly walks units in contig order. This
    is the fix for "thread-local concatenation alone is insufficient":
    concatenating per-thread buffers would make the order depend on the
    schedule.
  - Memory per unit is bounded by the unit's records plus one chunk of slack.
- **Window-based packing across scan boundaries.**
  - Windows carry rolling state but **no bases**: the next read starts exactly
    at the old `end` (:196-207).
  - RC packing reaches 28 bases left of `j`; forward reaches 40 bases right of
    `j`.
  - Fix: prepend a 39-base carry (the last `KMER-1` bases of the previous
    window) when advancing. Size `SEQ_buf` as `SCAN_MAX + 8 + 40` (+1 sentinel).
  - A seed whose forward k-mer crosses the not-yet-read right edge falls back
    to `Get_Contig_Piece`. That is at most a few seeds per window, and they are
    counted.
  - Packing reads the same NUMERIC bases the scan reads, so there is no
    encoding mismatch. Direct bit-copy from `.bps` 2-bit bytes is *not*
    planned: `.bps` puts the first base in the low bits (`Compress_Read`,
    `lib/gene_core.c:336`) while records pack it first in the high bits, so it
    would need a reversal table for little gain.
- **Thread-local buffers and peak memory.** Today's peak per rank is
  `local` (up to 2N with doubling slack) + `send_buf` (N) + `recv` (N), i.e.
  ~3–4N, where N is records × 16 B. Each step below states its peak.
- **EXTERNAL-mode validation.**
  - Unit tests use `fmemopen`-backed synthetic GDBs. That is EXTERNAL mode,
    but in memory, with no `seqpath`.
  - Every step also runs on real on-disk `.1gdb`/`.bps` files: the small
    fixture `tests/fixtures/generated/small6.1gdb` in a new unit/integration
    test, and yeast7 through `verify-phase1/2` and e2e.

## 5. Roadmap (revised order and scope)

Each step is a separate branch and commit, gated by the full suite: `make test`,
`make test-asan`, `test-mpi`, `verify-phase1/2/3-stage1/3-stage2/3-stage2-lcp`,
small e2e, and yeast7 PAF `cmp` at N = 1/2/4/8. Measurements always use the
Job D configuration: 1 node, `-c 16`, `--cpu-bind=cores`, `OMP_NUM_THREADS=1`,
3 reps, interleaved.

### O1-0 — Instrumentation and the tie-order question (new, first)

- **Code**
  - Per-rank `extract_s`, record count and `Get_Contig_Piece` call/byte
    counters, printed per rank. Today only the `MPI_MAX` is reported.
  - Split `sort_s` into `exchange_s` (histogram + Alltoall/Alltoallv) and
    `local_sort_s`.
  - Behind an env var only: `PGA_PRESORT_SHUFFLE=seed` permutes `recv` before
    `sort_records`.
- **Purpose**
  - Confirm or refute the I/O hypothesis before building on it.
  - Answer invariant 5: if yeast7 PAF is identical under 3 shuffle seeds at
    N = 1 and 4, then destination-internal order is free and O1c may use
    bucket-major order. Otherwise O1c must reproduce today's order exactly.
- **Benefit:** none (measurement only).
- **Risk:** timer placement only; default paths are byte-identical.

### O1a — Pack from the scan window (smallest independently testable change)

- **Code:** `src/extract.c` (39-base carry, window-relative packing, fallback
  counter) and new window-based variants of `pack_forward_kmer`/`pack_rc_kmer`
  in `src/pack.c`. The old functions stay for the fallback and as the test
  oracle.
- **HySortK reference:** encoding from resident sequence (`copy_bits`). The
  technique itself is generic.
- **Expected benefit:** removes 1–2 seek+read calls per record (~128M at N=1
  on yeast7). If O1-0 confirms the hypothesis, `extract_s` should drop
  substantially, and the N=8 anomaly (shared-file contention) should shrink.
- **Risks:** off-by-one at the carry and window edges, and the RC left reach at
  the window start.
- **Tests**
  - New unit test: window packing vs `pack_*_kmer` for every seed of E-06 and
    E-11.
  - `test_extract_window` (`SCAN_MAX=64`, which forces many window seams) plus
    a new seam-straddling case.
  - Byte-identical `extract_to_records` buffer vs the old path on
    `small6.1gdb` (on disk).
  - `verify-phase1` (byte-identical to `gold.tuples`) and all gates.
- **Measure:** per-rank `extract_s`, fallback count, syscall counters, RSS.

### O1b — Parallel scan over work units, deterministic assembly

- **Code**
  - `extract_to_records_mt` with a pthread pool. `PGA_EXTRACT_THREADS`
    defaults to 1 = today's path.
  - Per-thread GDB reader (§4); `Comp[]` initialised once.
  - Work units are whole contigs, claimed largest first via an atomic counter.
  - Each unit writes its own segment; assembly concatenates segments in contig
    order into `local`.
- **HySortK reference:** none beyond "threads write private slots"
  (`ScatteredSupermers[tid]`). Our determinism scheme differs.
- **Expected benefit:** extraction speedup bounded by the largest contig per
  rank and by cores per rank. With 8 ranks × `-c 16` the cores are otherwise
  idle in this stage.
- **Peak memory:** today's plus one copy during assembly (segments + `local`).
  O1c removes this.
- **Risks:** races (run TSan), a nondeterministic schedule leaking into the
  output order, and error paths inside threads.
- **Tests**
  - Byte-identical pre-sort buffer for T ∈ {1,2,3,4,8} vs T=1, on synthetic GDBs
    and on `small6.1gdb`.
  - New `make test-tsan` target.
  - yeast7 `pga-extract` dump identical across T; all gates at T ∈ {1,4}.
- **Measure:** extraction strong scaling in T at fixed N, and ranks × threads
  at a fixed core count (8×1 vs 8×4 vs 2×16).

### O1c — Destination-aware staging (the HySortK adaptation)

Replaces `local` → `bucket_local_counts` (a separate full pass) →
`organize_send_buffer` (count pass + per-record `memcpy`) with:

1. Each work unit emits into **per-bucket chunked segments**: 1024 routing
   buckets, playing the role of HySortK's tasks.
2. The bucket histogram is accumulated during emission, which removes the
   separate pass.
3. `MPI_Allreduce` → the unchanged `compute_destination_map` (contiguous
   ranges).
4. Exact per-destination offsets come from the per-(unit, bucket) counts.
   The send buffer is filled once, deterministically. Chunks are freed as they
   are consumed.

Candidate designs:

| Design | Scans | Copies into send_buf | Peak (records) | Sync | Order |
|---|---|---|---|---|---|
| Today | 1 + histogram pass + count pass | 1 | ~3–4N | 1 Allreduce | extraction order per dest |
| A. Two-pass: count-only scan → Allreduce → fill scan writing straight into send slots | 2 | 0 (written in place) | 2N | same | can reproduce today's order exactly |
| **B. Unit × bucket staging → Allreduce → deterministic pack** (recommended) | 1 | 1 | ~2N + chunk slack | same | (dest, bucket, unit, scan) by default; today's order needs a unit-major walk |
| C. Today's `local` + fused histogram only | 1 | 1 | ~3–4N | same | unchanged |

- **Recommendation: B.**
  - It removes two full passes and the 2N doubling slack, and it is the
    faithful HySortK structure.
  - It needs a decision on pack order. If O1-0 shows order is free, pack in
    (bucket, unit) order. If not, walk (unit, then the buckets belonging to
    dest), which reproduces today's order at the same cost because per-unit
    bucket offsets are known.
- **Fallback: A**, if memory is the binding constraint at larger datasets and
  O1a/O1b made the scan cheap enough to run twice. Decide from measured scan
  cost, not by assumption.
- **Do not assume fusion is faster:**
  - The removed passes are streaming passes over N × 16 B, cheap relative to
    the current extraction cost.
  - Per-bucket chunks add scattered writes across 1024 streams per unit; mind
    TLB and cache effects, so keep chunks ≥ 64 KB and a per-thread chunk pool.
  - The clearest win is **memory** (~4N → ~2N). Runtime gains must be
    measured.
- **Code:** a new `src/stage_buf.{c,h}`; `run_stage2` gains an entry point
  that takes staged input (the current one is kept for `test_stage2`);
  `pga-mpi-merge.c` calls the new path.
- **Risks:** offset arithmetic, `INT_MAX` guards (keep the existing ones), and
  order equivalence.
- **Tests**
  - Per-destination send byte streams identical to `organize_send_buffer` on
    the same input (or identical as a multiset, if O1-0 frees the order), for
    S-09 inputs and for `small6.1gdb` at simulated world sizes 1–8.
  - `test-mpi` modes uniform/skew/allempty/seam; all gates.
- **Measure:** `extract_s + exchange_s`, peak RSS (`/usr/bin/time -v` per
  rank), and pass counts.

### O1d — Batched nonblocking exchange (deferred; belongs with O2a)

- **Expected benefit at the current scale is small.** The whole exchange plus
  local sort is 1–4 s against 20–49 s of extraction. Overlap can hide at most
  the exchange time.
- **MPI progress.** cray-mpich makes no asynchronous progress by default.
  HySortK's overlap works because packing happens between `Ialltoall` and
  `Wait`. For us, overlapping with *extraction* would require interleaving
  batches with the scan, i.e. sends driven by fill level. That is a bigger
  restructuring (and threads plus MPI need `MPI_THREAD_FUNNELED`).
- **What it would genuinely solve:**
  - The `int` count/displacement ceiling (O2a): fixed 80 KB slots never
    approach `INT_MAX`.
  - Memory: streaming from staging removes the full `send_buf`, giving
    ~N staging + N recv.
- **Cost:** a variable number of rounds, padding of empty slots, done-flag
  rounds, and determinism of the receive order. Receive order must be source
  order, then the sender's order; that is achievable because slots are
  per-source.
- **Recommendation:** after O1c, measure `exchange_s` at multi-node scale (2–4
  nodes). Implement only if exchange is a significant fraction, or as the
  O2a large-count fix. Keep it out of O1.

## 6. Supermer compression: not applicable

- Supermers pay off because *adjacent k-mers share a minimizer and therefore a
  destination*, so a run of ~10–200 k-mers is sent as one 2-bit string.
- Our destination is the 5-base prefix of each seed's own k-mer. Adjacent seeds
  are uncorrelated, so the expected run length at P ranks is ~P/(P−1), i.e.
  ≈1. Seeds are also dense here: ≈0.77 records per base on yeast7.
- Making supermers apply would require minimizer routing, which breaks
  prefix-range ownership (§1). Rejected.
- **Weak analogue, optional:** the LCP and mask bytes are always 0 before
  sorting, so a sender could ship 14 of 16 bytes. The receiver re-expands
  directly into the 16-byte layout during the receive copy, which is HySortK's
  "expand in the receive loop" idea.
  - That is 12.5% less exchange traffic, but it adds a copy and the exchange is
    1–4 s.
  - Only worth it inside O1d. Not planned otherwise.

## 7. Recommendation summary

1. Finish the gate (C2, D, `baseline.md`).
2. O1-0, the instrumentation and tie-order experiment, is small and decides
   O1c's design.
3. O1a is the smallest independently testable optimisation and likely the
   largest single win. It is generic, and should be labelled as such.
4. O1b, deterministic per-unit assembly.
5. O1c is the HySortK contribution: destination-aware staging with a fused
   histogram, ~2N peak memory. Report memory and runtime separately.
6. O1d is deferred to O2a and multi-node measurements.
