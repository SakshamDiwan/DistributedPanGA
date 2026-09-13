# Architecture

Source commit inspected: `ce64406b45908451bd161466ebab0646cc84b92d`.
Symbols below were read from the executable code, not from comments; where a
comment contradicts the code, the code is described.

## Phase map and entrypoints

| Phase | Binary | Driver | Adds |
|---|---|---|---|
| 1 | `pga-extract` | `src/pga-extract.c` | serial extract → text tuples |
| 2 | `pga-sort` | `src/pga-sort.c` | `sort_records` (msd_sort) |
| 3 stage 1 | `pga-mpi-stage1` | `src/pga-mpi-stage1.c` | gather-to-rank-0 baseline |
| 3 stage 2 | `pga-mpi` | `src/pga-mpi.c` | `run_stage2` distributed sort |
| 4 | `pga-merge` | `src/pga-merge.c` | single-rank adapter + FastGA merge/align |
| 4b | `pga-mpi-merge` | `src/pga-mpi-merge.c` | distributed merge + alignment (**final**) |

`pga-mpi` is the distributed-sorting diagnostic; it stops after the sort and
dumps tuples. `pga-mpi-merge` is the end-to-end executable and takes exactly
one argument, `<gdb-stem>` — no thread or frequency flags.

## Stage 1: extraction

`extract_to_records(gdb, num_workers, worker_id, buf)` → `src/extract.c:314`.

- `compute_work_split` (`src/work_split.c:3`) splits **whole contigs** across
  workers by cumulative base count: worker `w` owns contig indices
  `[split[w], split[w+1])`. No within-contig splitting, so a single dominant
  contig cannot be shared (this is proposal O2).
- `scan_contig` (`src/extract.c:65`) mirrors FastGA's `scan_thread`. It walks
  `SMER=8`-mer minimizers via `TMap`/`Comp` (`src/tables.c`) to find
  `TMER=12` syncmers, in `SCAN_MAX = 10_000_000`-base windows.
- Selection branches: `mz < min4` (new right-end), `pos4 == i - SOFF`
  (left-end rescan), `mz > min4` → `continue`. The `mz == min4` tie **falls
  through and emits**. Preserve all four cases.
- Emission is via an `EmitFn` callback, so Phase 1 (`Tuple`) and Phase 2+
  (`RecordBuffer`) share one scanner.
- Per emission point `j`: forward k-mer emitted when `j <= len - KMER`
  (`pack_forward_kmer`, position stored = `j`, strand 0); reverse-complement
  emitted when `j >= KMER - TMER` (`pack_rc_kmer`, position stored = `j+TMER`,
  strand 1). Both may fire at the same `j`.
- Extraction is **serial within each rank** — the only parallelism is across
  ranks.

### Packing (`src/pack.c`)

- Forward covers bases `[j, j+40)`; each output byte packs 4 bases MSB-first,
  `(b0<<6)|(b1<<4)|(b2<<2)|b3`.
- Reverse-complement reads forward `[j-28, j+12)` and, for output byte `b`,
  packs local indices `[36-4b, 39-4b]` then maps through `Comp[]`.
- Both use a 65-byte stack buffer passed as `raw+1` because
  `Get_Contig_Piece` writes a sentinel at `buffer[-1]`. Out-of-range returns
  `-1` and the emission is skipped.

## Record format (`src/record.h`)

```
offset 0                      LCP byte        (semantic, [0,40])
offset 1 .. KBYTES            packed 40-mer   (10 bytes)
offset 1+KBYTES               mask byte       (always 0; masking unsupported)
offset POS_OFFSET             position        (post_bytes, little-endian)
offset POS_OFFSET+post_bytes  contig+strand   (cont_bytes, little-endian)
```

`record_size = 1 + KBYTES + 1 + post_bytes + cont_bytes`.
`compute_record_sizing` derives `post_bytes` from `gdb->maxctg` and
`cont_bytes` from `2 * gdb->ncontig` (the factor 2 reserves the strand bit).
Strand lives in the **high bit of the last contig byte**
(`1 << (8*cont_bytes - 1)`). For yeast7: `post_bytes=3, cont_bytes=1,
record_size=16`.

Buffers are allocated as `capacity * record_size + 1`; the trailing byte is
`msd_sort`'s sentinel write at `array[asize]`. Growth doubles capacity and
preserves the `+1`.

## Stage 2: distributed sort (`src/dsort.c`)

`run_stage2(local, recv, comm)` → `src/dsort.c:312`:

1. `bucket_local_counts` — `record_bucket` = first 5 bases =
   `kmer[0]<<2 | kmer[1]>>6`, giving `NUM_BUCK=1024` buckets.
2. `MPI_Allreduce(MPI_LONG_LONG, MPI_SUM)` over the 1024 counts → identical
   global histogram on every rank.
3. `compute_destination_map` — prefix sums, then FastGA's `distribute()`
   splitter logic; fills `select[bucket] -> rank` and `ksplit[]`. Deterministic
   given identical input, so no communication is needed to agree on ownership.
4. `organize_send_buffer` — two passes (count, then place) into one contiguous
   destination-ordered buffer. Guards both per-rank and cumulative byte counts
   against `INT_MAX` and aborts.
5. `MPI_Alltoall` on the counts, then `MPI_Alltoallv` with `MPI_BYTE`. Receive
   displacements are guarded against `INT_MAX`.
6. `sort_records(recv)` — local radix sort (below).
7. `recalc_all_lcps(recv)` under `-DLCPs` — overwrites byte 0 buffer-wide.
8. `verify_per_rank_sorted` → abort code 2 on violation.
9. `verify_cross_rank_boundaries` — `MPI_Sendrecv` of `KBYTES` edge k-mers with
   `MPI_PROC_NULL` at the ends, `MPI_Allreduce` of the local verdict → abort
   code 3.
10. Seam fixup: for `my_rank > 0 && recv->count > 0`, byte 0 of record 0 is
    rewritten as `simple_kmer_lcp_bases(prev_last_kmer, my_first)`.

MPI calls used in the sort path: `Allreduce`, `Alltoall`, `Alltoallv`,
`Sendrecv`, `Abort`. Stage 1's `mpi_gather_records` uses `Gather`/`Gatherv`.

### Local sort (`src/sort.c`)

`sort_records` pre-buckets records by k-mer byte 0 into a temporary buffer,
builds `part[256]` byte lengths, then calls
`msd_sort(data, count, rsize, ksize = 1 + KBYTES, part, 0, 256, nthreads)`
with **`nthreads = 4` hardcoded**. `OMP_NUM_THREADS` has no effect. `msd_sort`
starts at digit 1, so the LCP byte at offset 0 is skipped as a key byte.

### LCP semantics

`simple_kmer_lcp_bases(a, b, kbytes)` returns matching bases in `[0, 40]`.
`recalc_all_lcps` writes 0 for record 0 and the pairwise value elsewhere
(`-DLCPs`), or a 0/1 boundary marker without it. This deliberately replaces
`msd_sort`'s byte 0, whose encoding mixes intra- and inter-bucket formulas and
disagrees with bases-in-common at `part[]` boundaries. GIXshow's `lcp` column
is a verbatim view of that raw byte and is **not** a valid oracle.

Callers matter: `pga-sort` calls `sort_records` **only**, so its byte 0 is raw
`msd_sort` output. Only `run_stage2` produces semantic LCPs.

## Stage 3: FastK/FastGA adapter (`src/kmer_adapter.c`)

`build_kmer_adapter(buf, T, P)` builds a `KmerStreamAdapter` whose field order
mirrors `Kmer_Stream` (`lib/libfastk.h`) and a `PostListAdapter` mirroring
`Post_List` (`src/fastga_pipeline.h`). Both are **cast** to the FastGA types at
the call site, so field order and widths are load-bearing.

- `ibyte = 3` → prefix is k-mer bytes 0..2, held implicitly in `cpre`;
  `hbyte = kbyte - ibyte = 7` suffix bytes are stored per entry.
- `pbyte = hbyte + 1 + 1 + post_bytes + cont_bytes` (= 13 for yeast7);
  `tbyte = ibyte + pbyte` (= 16) is the width `Current_Entry` emits.
- Table layout per entry: 7 suffix bytes, `0` (mask), the record's LCP byte,
  then position+contig bytes copied verbatim.
- `index[]` has `ixlen = 2^24` `int64` entries ≈ **128 MiB per rank**, allocated
  before records and cache. Semantics are cumulative **`prefix <= p`**.
- `P->maxp` = longest run of one 3-byte prefix `+ 1`; it sizes the merge cache
  (`NTHREADS*(maxp+1)*KBYTE`).
- Stream ops are redirected by macros at `src/fastga_pipeline.c:60-66`
  (`First_Kmer_Entry`, `Next_Kmer_Entry`, `GoTo_Kmer_Index`,
  `Clone_Kmer_Stream`, `Free_Kmer_Stream`, `Current_Entry`).
- Ownership: `free_kmer_adapter` frees `table` and `index`. Clones share both;
  `adapter_Free_Kmer_Stream` frees only a clone's wrapper. End of stream is
  `csuf == NULL` with `cpre == ixlen`.
- `adapter_GoTo_Kmer_Index` is an O(ixlen) scan, tolerable only because
  `self_adaptamer_merge` calls it `NTHREADS-1` times; with `NTHREADS=1` it is
  never called from the split loop.

## Stage 4: contig-pair assignment (`src/contig_assignment.c`)

`contig_assignment_build(gdb, world_size, out)`:

- Enumerates all `N*(N+1)/2` pairs `(i,j), i <= j`, weight `clen[i]*clen[j]`.
- `qsort` with `pair_cmp_weight_desc`: weight descending, then `i` then `j`
  ascending — deterministic across libcs.
- Least-loaded-rank assignment via a min-heap keyed `(load, rank)`; ties go to
  the lower rank id. Every rank computes the same map with no communication,
  relying on replicated GDB metadata.
- Storage: `int8_t rank_for_pair[N*N]` with only the upper triangle written and
  `-1` elsewhere → quadratic memory, and `world_size > 127` is rejected.
- `pga-mpi-merge.c:406` then builds a **second** `NCONTS*NCONTS` symmetric map
  (`RankForPair`) so the merge inner loop can index `[i*NCONTS + j]` without a
  min/max swap. Two quadratic maps are live simultaneously.

Worked example (independent of the code): two equal-length contigs on two
ranks. All three pairs tie on weight, so `(i,j)` order is
`(0,0), (0,1), (1,1)` and the heap yields ranks **0, 1, 0**.

## Stage 5: merge, seed exchange, alignment

`pga-mpi-merge.c:main` sequence, with the stage timers it reports:

- **extract** — `record_buffer_init` + `extract_to_records`.
- **sort** — `run_stage2`.
- **merge** — `build_kmer_adapter`, FastGA extern setup (`IBYTE`, `ICONT`,
  `IPOST`, `ISIGN`, `KBYTE`, `CBYTE`, `LBYTE`, `PAYOFF`, `ESHIFT`, `NCONTS`,
  `AMXPOS`, `BMXPOS`, `MAXDAG`, `DBYTE`, `Perm1`, `IDBsplit`, `Select`),
  `contig_assignment_build`, `RankForPair`, per-rank scratch dir, `IOBuffer`
  allocation, `self_adaptamer_merge`, both `redistribute_seeds` calls, and both
  `recompute_buck` calls.
- **alignments** — `pair_sort_search` **and** the `ALNtoPAF` conversion.

Driver-set externs of note: `NTHREADS = 1`, `NPARTS = 1`, `SELF = 1`,
`FREQ = 10`, `CHAIN_BREAK = 2000`, `CHAIN_MIN = 170`, `ALIGN_MIN = 100`,
`ALIGN_RATE = 0.3`. `pga-merge.c` sets the same `NTHREADS=1`/`FREQ=10`.

### Seed emission and routing

`new_self_merge_thread` (`src/fastga_pipeline.c:231`) walks adaptamer groups
and, for each seed pair, computes the destination:

```c
idest = Select[icont];                       /* single-rank fallback */
if (RankForPair != NULL)                     /* distributed mode */
    idest = RankForPair[icont * (int64) NCONTS + jcont];
ou = (isign == jsign) ? nunit + idest : cunit + idest;
```

`jcont` is masked with `((1 << (8*JCONT - 1)) - 1)` to strip the strand bit.
Each seed is `1 + IBYTE + JBYTE` bytes: `plen`, then `pay1`, then `p`.
Frequency cutoff is `kfreq = FREQ * kbyte`. Buffers flush to
`$SORT_PATH/$PAIR_NAME.<k>.{N,C}` when `btop >= bend`.

### Seed exchange (`redistribute_seeds`, `src/pga-mpi-merge.c:132`)

Per seed type (`N` = same-strand, `C` = cross-strand):

- `world_size == 1` shortcut: keeps `units[0].file` and returns its size
  without any MPI call. This dodges `INT_MAX` on one rank only — the multi-rank
  path still has count limits.
- Reads each per-destination file into one send buffer, `close`s and `unlink`s
  the originals, `MPI_Alltoall`s the counts, then `MPI_Alltoallv`s `MPI_BYTE`.
- Writes the received bytes to `$SORT_PATH/$PAIR_NAME.combined.{N,C}` and
  installs that fd as `units[0].file`, then sets `NUNITS = 1`.
- Guards: per-destination send size vs `INT_MAX`, and cumulative **receive**
  total vs `INT_MAX`. There is **no** guard on the cumulative *send* total
  before `send_displs[r] = (int) cum` — see `state.md` follow-ups.

`recompute_buck` then rescans the combined file to rebuild per-A-contig seed
counts (`unit->buck[icont]`), because redistribution invalidates the counts the
merge wrote. It reads the contig field at offset `1 + IPOST` and strips the
strand bit.

### Alignment and output

`pair_sort_search(gdb, gdb)` sorts and chains seeds, writing
`$SORT_PATH/$ALGN_UNIQ.<tid>.las` / `$ALGN_PAIR.<tid>.las` scratch, and
`la_merge` produces `$SORT_PATH/$ALGN_UNIQ.1aln`
(`src/fastga_pipeline.c:1369`). The driver then shells out:

```c
snprintf(cmd, ..., "ALNtoPAF %s > out.rank%d.paf 2>/dev/null", aln_path, my_rank);
```

## Disk interfaces and lifetimes

| Path | Written by | Scope |
|---|---|---|
| `<gdb-stem>.1gdb`, `.gix` | FastGA tools (external) | read-only input, must be visible to all ranks |
| `/tmp/pga-mpi-merge-r<rank>-XXXXXX` | `make_temp_dir_for_rank` | **node-local**, hardcoded `/tmp`; `TMPDIR` is ignored |
| `$SORT_PATH/_pair.<pid>.<k>.{N,C}` | merge threads | per-rank scratch |
| `$SORT_PATH/_pair.<pid>.combined.{N,C}` | `redistribute_seeds` | per-rank scratch |
| `$SORT_PATH/_uniq.<pid>.1aln` | `la_merge` | per-rank scratch |
| `out.rank<N>.paf` | `ALNtoPAF` via `system()` | **current working directory** |
| `pga-mpi-merge.timing.csv` | rank 0 | current working directory, appended |
| `<out-stem>.<rank>.tuples` | `pga-mpi` | current working directory / given path |

Scratch is removed by `system("rm -rf $SORT_PATH")` at exit. Because outputs
are CWD-relative and scratch is `/tmp`-local, multi-node runs need a unique
**shared** CWD and give no shared visibility of scratch.

## Timing CSV

Columns: `world_size,extract_s,sort_s,merge_s,alignments_s,total_s`. Each stage
value is `MPI_MAX` across ranks (`MPI_Reduce`). **`total_s` is the sum of those
four maxima**, not an independently measured end-to-end wall time; different
ranks can be the max in different stages, so it overestimates. Measure real
wall time separately. The header is written only when the file is empty, and
rows are appended.

## Text dump formats (`pga-mpi`, `pga-mpi-stage1`, `pga-extract`)

- Default, 4 columns: `<contig> <position> <strand> <hex_kmer>` — k-mer is `$4`.
- `PGA_DUMP_LCP=1` (`pga-mpi` only), 5 columns:
  `<contig> <position> <strand> <lcp> <hex_kmer>` — LCP is `$4`, k-mer is `$5`.

Scripts must key on the right field for the mode they run in.
