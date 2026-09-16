# Test catalog

Fast checks on the current pipeline so the O1–O5 work can't regress it
silently. The pre-existing safety net is byte-identical comparison against
FastGA, which needs gold fixtures, the FastGA tools and an allocation — none of
which exist in this checkout.

Tests encode **intended contracts**. Expected values come from headers, the
record-layout spec and the algorithm definition — never captured from a run.
Where the code deviates from its own contract that becomes a regression test
asserting the *correct* answer (`R-nn`), marked `EXPECT_FAIL_UNTIL` on the one
failing assertion. A known bug never becomes a passing expectation.

## Status

Current results live in `tests/RESULTS.md`. Summary: **42,471 checks,
0 failures, 5 xfail**, `make test` ≈ 25 s, no allocation or fixtures needed.

| Group | IDs | Binary | State |
|---|---|---|---|
| Extract | E-01…E-11, E-13 | `test_extract`, `test_extract_window`, `test_record_split` | done, passing |
| Work split | E-12 | `test_record_split` | done, passing |
| Sort / LCP | S-01…S-05b | `test_lcp_sort` | done, passing |
| Splitter / exchange | S-06…S-09 | `test_dsort` | done, passing |
| Pair assignment | A-01…A-05 | `test_assign` | done, passing |
| Adapter | K-01…K-05 | `test_adapter` | done, passing |
| Real MPI | M-01…M-06 | `test_stage2` | written; needs an allocation (world=1 smoke-tested) |
| Regression | R-02 | `test_exchange_plan` | **fixed**, markers removed |
| Regressions | R-03, R-04 | — | done, **xfail** (4 + 1 assertions) |
| Regression | R-01 | `test_stage2` | written; needs an allocation |
| End-to-end | Z-01, Z-02 | — | blocked: no fixtures, no FastGA tools |

## Scope

Latest stage of each phase only: `extract_to_records` + `pack.c`; `run_stage2`;
`build_kmer_adapter` and pair assignment; one end-to-end comparison.

**Not tested:** Slurm/Perlmutter mechanics; shell-script plumbing; vendored
FastGA in `lib/` (wave aligner, chaining, `msd_sort` internals — preserving that
boundary is what makes byte-identical validation mean anything); superseded
phases (`pga-extract`, `pga-sort`, `pga-mpi-stage1`, `extract_to_tuples`,
`mpi_gather_records`); merge internals (`self_adaptamer_merge`,
`pair_sort_search`, `la_merge` — covered by Z-02; `redistribute_seeds` and
`recompute_buck` are `static` in a file with `main()`, so unit-testing them
means moving them to `src/seed_exchange.c`, which belongs to O4); O5/GPU.

## Conventions

**MPI environment.** The default Perlmutter environment loads
`craype-accel-nvidia80`, so cray-mpich aborts with "GPU_SUPPORT_ENABLED is
requested, but GTL library is not linked". This is a CPU-only project:
`tests/run_mpi.sh` exports `MPICH_GPU_SUPPORT_ENABLED=0` (`module load cpu`
also works).

**Flags.** Test targets use `$(CFLAGS)`, so `-DLCPs` matches the `LCPS=1`
default. A test built without it silently exercises `recalc_all_lcps`'s 0/1
branch instead of the semantic `[0,40]` one — this was hit during
implementation and cost a false failure.

**Harness.** `tests/support/check.h|.c`: `CHECK`, `CHECK_MSG`, `CHECK_EQ_I`,
`EXPECT_FAIL_UNTIL(id, expr)`, `RUN(fn)`. No framework, no registration magic;
each binary has an explicit `main()`. XPASS *fails* the suite ("bug fixed,
remove the marker") so a marker cannot outlive its fix. Crashes, aborts and
timeouts are never acceptable expected failures.

**Linking.** One binary per concern, each linking only what it needs —
`test_record_split` links just `record.c` and `work_split.c`; only the
extraction binaries pull in `GDB.c`. Never `align.c`, `alncode.c`, `libfastk.c`
or `RSDsort.c`.

**Synthetic GDB** (`tests/support/synth_gdb.c`). `Get_Contig_Piece`
(`lib/GDB.c:1831`) reads only `seqs`, `seqstate`, `ncontig` and
`contigs[i].{boff,clen}`, so a GDB is built in memory with
`seqstate = EXTERNAL` over an `fmemopen`'d 2-bit blob — byte-for-byte the path
production takes. Verified to round-trip exactly, including unaligned `beg` and
both sentinels.

- Contig `i` starts at byte `boff`, occupies `(clen+3)>>2` bytes, byte-aligned.
  `Compress_Read` builds the blob so `Uncompress_Read` inverts it exactly.
- **Never** use `seqstate = NUMERIC`/`COMPRESSED`: those in-memory branches
  (`lib/GDB.c:1868-1895`) compute the source as `boff + beg/4 + beg%4`, wrong
  for `beg >= 4` — exactly how `pack_forward_kmer` is called. Production never
  reaches them.
- Always pass `buf + 1`; `buffer[-1] = 4` is a real write.
- Never `Close_GDB` a synthetic GDB — it would free fields we never allocated.

**Contig lengths.** Forward needs `j <= len-40`, RC needs `j >= 28`, with
`j <= len-12`. So `len >= 40` permits both orientations at *different* offsets;
`len >= 68` is needed for both at the *same* offset. Max stored position is
`len` (RC stores `j+12` with `j` up to `len-12`).

## Prerequisites

| ID | Change | For | State |
|---|---|---|---|
| PR-1 | `SCAN_MAX` wrapped in `#ifndef` (`src/extract.c`) | E-11 | **done** — 5 lines, zero runtime effect |
| PR-2 | Move the count/displacement arithmetic out of `redistribute_seeds` into a non-static, MPI-free `plan_byte_exchange()` in `src/seed_exchange.c` | R-02 | **done** |
| PR-3 | `.gitignore`: `tests/bin/` + `tests/fixtures/generated/`, `!tests/**/*.md`; drop the orphan `test_contig_assignment` line | all | **done** |

```c
/* PR-2. 0 on success; <0 if any per-slot count, or the cumulative total,
   would overflow the int used by MPI_Alltoallv counts/displacements. */
int plan_byte_exchange(const int64_t *sizes, int n, int *counts, int *displs);
```

`organize_send_buffer` (`src/dsort.c:189-202`) already guards both per-slot and
cumulative totals; it is the model. Only the `pga-mpi-merge.c` path is missing
the cumulative *send* guard.

---

## P1 · Extract

`scan_contig`, `record_emit_cb` and `EmitFn` are private to `src/extract.c`, so
every test drives public `extract_to_records(gdb, 1, 0, &buf)` and decodes with
the independent decoder in `tests/support/rec_build.c`.

| ID | Pre | Post |
|---|---|---|
| E-01 | `"ACGT"`×13 | `pack_forward_kmer(…,0,…)` → 0, every byte `0x1B` (A=0,C=1,G=2,T=3 → `00011011`b) |
| E-02 | 40 `A`s | `pack_rc_kmer(…,28,…)` → 0, every byte `0xFF`; independently `Comp[0x00] = 192\|48\|12\|3` |
| E-03 | 100 b random | forward and RC packing match the independent packers across many offsets |
| E-04 | 52 b | `j=12` → 0 (`[12,52)` fits exactly); `j=13` → `-1` |
| E-05 | 52 b | `j=27` (`begb=-1`) and `j=41` (`endb=53`) → `-1`; `j=28`, `j=40` → 0 |
| E-07 | 2 kb random | every forward record's position is a selected offset `j <= len-40`; every RC position is `j+12` with `j >= 28` |
| E-08 | 41 b | forward only at `j <= 1`, RC only at `j >= 28`; no offset yields both (1 < 28); full multiset matches oracle |
| E-09 | 67 b and 68 b, `j=28` selected | at 67 → RC at 40 only (28 > 27); at 68 → **both** forward at 28 and RC at 40 |
| E-10 | 128 contigs × 80 b | `cont_bytes == 1`; every decoded contig in range; contig 127 (max id) yields both strands with the id uncorrupted |
| E-12 | 8 equal / 1 giant+20 small / 1 contig / workers>contigs | `split[0]=0`, `split[nw]=ncontig`, non-decreasing, every contig owned exactly once, `post[nw]==seqtot`; empty ranges legitimate |
| E-13 | `maxctg=1.5M`, `ncontig=120` | `post=3`, `cont=1`, `record_size == 1+10+1+3+1 == 16` |

**E-06 · extraction matches an independent oracle — the O1 gate.**
`tests/support/ref_syncmer.c` implements the contract naively and separately,
sharing only `TMap[]` (constant data, part of the definition) and recomputing
the 4-base reverse complement itself:

```
mz(i) = min(fwd, rc)
fwd   = (TMap[b(i..i+3)] << 8) | TMap[b(i+4..i+7)]
rc    = (TMap[rc4(b(i+4..i+7))] << 8) | TMap[rc4(b(i..i+3))]
for j in [0, len-12]:
    h[t] = mz(j+t), t in 0..4        # the 12-mer's 5 constituent 8-mers
    m = min(h);  a = earliest argmin
    select j  iff  a == 0 or h[4] == m      # closed syncmer, ties select
        j <= len-40 -> (pos=j,    strand=0, kmer=fwd40(j))
        j >= 28     -> (pos=j+12, strand=1, kmer=rc40(j))
```

Pre: one 4 kb fixed-seed contig. Post: occurrence-multiset equality on
`(contig, position, strand, kmer)`. **Anti-vacuity:** the oracle counts the
right-end, left-end-rescan and tie branches and asserts each `> 0`. Observed
`right=762 left=801 tie=2` over 1564 selected offsets / 3106 records — the tie
path fires twice, so without that counter a fixture could silently never reach
it.

**E-11 · multi-window scan.** `test_extract_window` is built `-DSCAN_MAX=64`
(PR-1), so a 200-base contig forces ≥ 4 windows and exercises the re-window
pointer arithmetic at `src/extract.c:196-203` (`Get_Contig_Piece(...) - beg`).
Both binaries run the *same* fixtures (`WINDOW_FIXTURE_SEEDS`) against the
*same* oracle, so single-window/multi-window agreement follows transitively —
and a seam bug corrupting both paths identically still cannot pass, because the
oracle is independent of both.

---

## P2 · Sort and semantic LCP

| ID | Pre | Post |
|---|---|---|
| S-01 | hand table + 4000 random pairs | `simple_kmer_lcp_bases` matches hand values and the independent `ref_lcp_bases`: identical→40, differ at base 0→0, 3→3, 4→4, 7→7, 8→8, 39→39. Range `[0,40]`, not `msd_sort`'s byte 0, not GIXshow's `lcp` |
| S-02 | records `k1,k1,k2` differing first at base 17, LCP bytes pre-set to 99 | `recalc_all_lcps` → `0, 40, 17` |
| S-03 | `count == 0` | no write, no crash |
| S-04 | 5000 records from a 1000-k-mer pool | adjacent k-mers non-decreasing; occurrence multiset of `(kmer,pos,contig,strand)` unchanged, **byte 0 excluded** (`msd_sort` rewrites it). Anti-vacuity asserts duplicates exist — observed 4007 duplicate-k-mer adjacencies, 194 identical occurrences |
| S-05a | hand-built, `count == capacity == 64`, allocated exactly `count*rsize + 1` | `sort_records` clean under ASan. **Scope limit:** the test allocates the `+1`, so it constrains the *sorter* only |
| S-05b | buffer from production `record_buffer_init` + `extract_to_records` on 1.5 Mbp, forcing a grow | probes `data[capacity*record_size]` — exactly the byte production promises. Observed `count=1190552 capacity=2097152`. **Verified to fail:** removing `+1` from both allocators trips ASan at that line |

Not yet implemented, in `test_dsort` (mpicc — `dsort.h` includes `mpi.h`):

- **S-06** `record_bucket` over all 256 `kmer[0]` × 4 top-bit patterns of
  `kmer[1]` equals `(kmer[0]<<2)|(kmer[1]>>6)`; 1024 distinct, covering
  `[0,NUM_BUCK)`; low 6 bits of `kmer[1]` irrelevant.
- **S-07** `bucket_local_counts` matches a hand tally, others 0, sum
  `== buf->count`.
- **S-08 — the O2c gate.** Histograms: uniform / one bucket / all zeros / two
  buckets, at world 1, 3, 4, 1024. `select` **monotone non-decreasing** (the
  invariant global ordering depends on); every `select[b] ∈ [0,world_size)`;
  `ksplit[world_size] == NUM_BUCK`, non-decreasing, consistent with `select`;
  deterministic across calls (this is what licenses computing it redundantly per
  rank with no communication); uniform at world 4 → `select[b] == b/256`.
- **S-09** `organize_send_buffer`: `displs[0]==0`,
  `displs[r+1]==displs[r]+counts[r]`, `sum(counts)==count*rsize`, every record
  present exactly once byte-identical, and a zero-record destination gets
  `count==0` with a defined displacement.

---

## P3 · Pair assignment (not started)

`src/contig_assignment.c`. Expected values from the documented policy
(`contig_assignment.h:15-17`): weight descending, `(i,j)` ascending tie-break,
least-loaded rank, lower rank id on load ties.

- **A-01** 2 contigs len `L`, world 2 → `(0,0)→0, (0,1)→1, (1,1)→0`.
- **A-02** 3 contigs len `L`, world 2 → pairs in `(i,j)` order → `0,1,0,1,0,1`.
- **A-03** 6 distinct lengths, world 3 → every `i<=j` owned in range, every
  `i>j` still `-1`, `sum(load_per_rank)` equals total weight, owner symmetric.
- **A-04** two builds byte-identical.
- **A-05** world 127 succeeds; world 128 returns non-zero and leaks nothing
  (ASan) — the `int8_t` cap.

## P4 · Adapter (not started)

`src/kmer_adapter.c`. Each adapter allocates 128 MiB for `index`, so keep these
few.

- **K-01** `index[p]` equals the count of **input** records with prefix `<= p`,
  counted by a naive loop over the input, never through the adapter;
  `index[ixlen-1] == nels`; non-decreasing.
- **K-02** `ibyte=3, kbyte=10, hbyte=7, pbyte=13, tbyte=16, ixlen=1<<24`.
- **K-03** entry layout: 7 suffix bytes, `0`, the record's LCP byte, then
  `post+cont` bytes verbatim from `POS_OFFSET`.
- **K-04** traversal visits exactly `nels` in order; `Current_Entry`'s first 10
  bytes (3 rebuilt from `cpre` + 7 suffix) equal the original k-mer; end of
  stream is `csuf == NULL && cpre == ixlen`.
- **K-05** clone shares `table`/`index` pointers, advancing it doesn't move the
  parent, free clone then parent → no leak or double-free (ASan).

## P5 · Real MPI (not started)

`run_stage2`, via `tests/run_mpi.sh` at `RANK_COUNTS="1 2 3 4 8"` — **3 matters**,
since the splitter's rounding and empty-slice paths only appear at awkward rank
counts. Per-rank `RecordBuffer`s are built by hand from a rank-seeded PRNG, so
no GDB and no fixtures. Gather with plain `MPI_Gatherv`, not
`mpi_gather_records` (stage-1 code). Exit codes 2 and 3 are `run_stage2`'s own
aborts and count as **failures**, never expected failures.

- **M-01** occurrence multiset of `(kmer,pos,contig,strand)` out equals in,
  byte 0 excluded.
- **M-02** each slice internally non-decreasing; for adjacent **non-empty**
  ranks `last(r) <= first(r')`. Empty ranks are skipped, not compared against
  their sentinels.
- **M-03** concatenate in rank order on rank 0 and recompute every LCP with the
  **independent** `ref_lcp_bases` — not production `simple_kmer_lcp_bases`,
  which would be self-referential. Record 0 of the concatenation is 0.
- **M-04** a rank starting with 0 records, and one with exactly 1; assert on the
  post-exchange slice, since a rank starting empty may still receive.
- **M-05** repeatability: identical seeds twice → the three semantic invariants
  hold identically. **Not** byte-identical slices: records sharing a k-mer may
  be permuted among themselves, and stable ordering of equal keys is not a
  documented contract of `sort_records`/`msd_sort`.
- **M-06** same global record set partitioned across world 1, 2, 3, 4 → identical
  occurrence multiset, globally non-decreasing k-mer order, and matching LCP
  bytes. **Not** required: identical ordering of distinct occurrences sharing a
  k-mer, since different world sizes feed different subsets to each rank's sort.
  Consistency evidence, not an independent oracle.

---

## Regression tests

### R-01 · LCP seam past an empty intervening rank (not started)

`verify_cross_rank_boundaries` (`src/dsort.c:253`), fixup at `:430-439`.

**Fixture**, derived from the splitter arithmetic: world 4; every record
`kmer[0] = 0xC0`; **4 records per populated bucket, 8 total**, distinct within
each bucket (vary bytes 2..9). Mass in bucket **768** (`kmer[1]>>6 == 0`) and
bucket **770** (`kmer[1]>>6 == 2`), leaving 769 empty. Tracing
`compute_destination_map` gives `select[768]=0, select[769]=1, select[770]=2`.

**Pre-assertions, all before the seam check:** those three `select` values;
`recv.count == 4` on ranks 0 and 2 and `== 0` on rank 1; ranks 0 and 2 are the
nearest non-empty neighbours. If the splitter changes these fail as
"precondition not met" rather than letting the seam assertion pass for the wrong
reason.

**Marked assertion:** rank 2's record-0 LCP equals
`ref_lcp_bases(last_kmer_of_rank_0, first_kmer_of_rank_2)` = **4** — both share
`kmer[0] = 0xC0` (T,A,A,A) and diverge in `kmer[1]`'s top bits.

**Today:** an empty rank sends `my_last = 0x00…` (`src/dsort.c:279`) and the
fixup compares against that sentinel instead of the true predecessor further
back; the `MPI_Sendrecv` only ever looks one hop. `LCP(0x00…, 0xC0…) = 0`. The
boundary check passes vacuously, so nothing aborts and the wrong LCP flows into
seed matching. **Fix:** carry the last k-mer of the last *non-empty* rank
(`MPI_Exscan`-style, or skip empty predecessors).

### R-02 · Cumulative send displacement overflow — **FIXED**

**Pre:** `sizes = {1<<30, 1<<30, 1<<30}` — each slot under `INT_MAX`, only the
running total overflows. Plain integers, so the test allocates nothing.
**Marked:** `plan_byte_exchange` returns non-zero. **Companions that must pass
today**, guarding the extraction: `{100,200,300}` → 0 with
`counts={100,200,300}`, `displs={0,100,300}`; one slot of `INT_MAX+1` → non-zero.

**Was:** the per-destination check and the cumulative *receive* check existed,
but `send_displs[r] = (int) cum` cast a running `int64` with no cumulative
*send* guard. Observed before the fix: `displs = {0, 1073741824, -2147483648}` —
the third wrapped negative, so MPI would have read outside the send buffer. The
`world_size == 1` shortcut skips MPI entirely, hiding this on one rank only.

**Fixed** by the cumulative guard in `plan_byte_exchange`. The marker was
removed only after XPASS flagged it, which confirmed the test had real power.

### R-03 · Position field representability — **confirmed empirically**

RC records store `j + TMER` with `j` up to `len - TMER` (`src/extract.c:189`),
so positions span `[0, maxctg]` — `maxctg + 1` values. `post_bytes` must
satisfy `maxctg <= 256^post_bytes - 1`. Expectations derived from
representability, **not** from production `bytes_needed`.

| `maxctg` | intended | current | |
|---|---|---|---|
| 1 | 1 | **0** | wrong (unreachable: no records at len 1) |
| 255 | 1 | 1 | ok |
| **256** | **2** | **1** | **wrong** |
| 257 / 65535 | 2 | 2 | ok |
| **65536** | **3** | **2** | **wrong** |
| 65537 | 3 | 3 | ok |

Fails exactly at `maxctg = 256^k`, where `while (cum < max_val)` stops a byte
short.

**Extraction half (marked):** a 256-base contig with `j = 244` selected — the
only offset that can yield RC position 256, since `j <= 244` and `j+12 >= 256`.
The test asserts that precondition, then asserts the decoded RC position is 256.
**Observed:** it decodes as **0**, truncated by `post_bytes = 1`.

**Passing companion:** `cont_bytes` is correct at its boundaries (128 → 1 byte,
129 → 2); the `2 *` factor in `bytes_needed(2 * ncontig)` already reserves the
strand bit, so only `post_bytes` is affected.

**Fix:** `bytes_needed(gdb->maxctg + 1)` or a representability-based helper.
`record_size` is an internal binary width — `gold.tuples` is text and
`gold.paf` is alignment output, neither encodes `post_bytes` — so this is **not**
a reason to regenerate independent gold. Regenerate only artifacts that truly
depend on the record format, and investigate any changed position or alignment
as a real behavioural change rather than refreshing the reference.

### R-04 · Singleton LCP byte

`src/sort.h:12-13` says "overwrite byte 0 of **every** record"; the first record
has no predecessor, so 0 is the only consistent value.

**Pre:** `count == 1`, byte 0 pre-set to 37. **Marked:** byte 0 becomes 0.
**Observed:** left at 37 — `if (buf->count <= 1) return;` (`src/sort.c:40`)
bundles the empty case, correctly a no-op, with the singleton case, which is
not. Latent today; the empty-rank work in R-01 makes single-record slices
reachable. **Fix:**
`if (count == 0) return; data[0] = 0; if (count == 1) return;`

---

## P6 · End-to-end (blocked)

Blocked on `tests/fixtures/`, `~/FASTGA` and `ALNtoPAF`, none of which exist
here.

Every existing script hardcodes the yeast fixture path and a FastGA location
(`runtime_test_phase4b.sh` points at another user's home). So this tier needs a
small parameterized runner, `tests/run_e2e.sh`, taking `--gdb-stem`,
`--gold-dir`, `--fastga-bin` (env `PGA_GDB_STEM`, `PGA_GOLD_DIR`,
`PGA_FASTGA_BIN`) — or the same three variables added to the existing scripts
with current values as defaults, whichever is the smaller diff. Keep it narrow
and separately reviewable.

Prefer a **~50 kb synthetic multi-contig FASTA** over 84 Mbp yeast7 for speed:
generator and seed tracked in `tests/fixtures/small/`, derived GDB/GIX/gold in
`tests/fixtures/generated/` (ignored). It must yield a **meaningful number of
alignments** — build in deliberate repeats and near-repeats, and assert a
positive lower bound recorded at gold-generation time, or every assertion passes
vacuously. Keep yeast7 as the slow authoritative check.

- **Z-01** `verify-phase1`, `verify-phase2`, `verify-phase3-stage2`,
  `verify-phase3-stage2-lcp` all PASS. For stage 2 **both** tiers are required:
  Tier 2 (k-mer order over the unsorted rank-ordered concatenation) proves the
  sort really distributed; Tier 1 (multiset equality vs gold) proves nothing was
  lost. Neither implies the other. Run the LCP target last — its script starts
  with `make clean`.
- **Z-02** actually run `pga-mpi-merge` at world 1, 2, 4. The **aggregate** PAF
  (all `out.rank*.paf` concatenated and canonically sorted) equals `gold.paf`,
  produced independently by FastGA — that is the real assertion; cross-world-size
  agreement alone is only consistency evidence. An individual rank file **may
  legitimately be empty** (a rank can own pairs yielding no alignment), so
  assert on the aggregate being non-empty, not on each file. A zero exit status
  does not mean PAF was produced: `ALNtoPAF` failure is printed and ignored
  (`src/pga-mpi-merge.c:519`). Do not cite the historical 43,126-alignment
  figure as reproduced — it belongs to a specific input, parameters and FastGA
  revision.

---

## Running

```bash
make test            # all unit binaries: no allocation, no fixtures
make test-asan       # ASan/UBSan, includes the slow production-growth test
make test-extract    # one concern
make test-clean

salloc --account=m4341 --constraint=cpu --qos=interactive --nodes=1 --time=00:30:00
make test-mpi                    # once P5 lands
RANK_COUNTS="1 3" make test-mpi   # awkward rank counts alone
```

**Acceptance:** `make test` green in seconds on a clean login node with no
fixtures and no FastGA tools; no unexplained failures — every reported failure
is either an `EXPECT_FAIL_UNTIL` marker naming a regression ID or a real
regression. Four is the number of known defect *categories*, not a target: the
marker count is whatever is currently unfixed and drops to zero as fixes land.
Pre-existing `make verify-phase*` behaviour unchanged.

**Suite self-checks**, so it cannot pass vacuously: E-06's branch counters prove
the tie and rescan paths were reached; S-04 asserts duplicates exist; R-01 and
R-03 assert their preconditions; XPASS fails the suite; S-05b was verified by
deliberately removing the `+1` and confirming ASan catches it.
