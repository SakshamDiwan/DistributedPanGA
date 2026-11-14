# O2 explained: large-count exchanges (`chunked_alltoallv`)

A walkthrough of the code behind report §4.2, for the author. Companion to
`o1-final-report.md` (O1) and `pr2-results.md` (O2 measurements). Source:
`src/chunked_exchange.h`, `src/chunked_exchange.c`, and the two callers in
`src/dsort.c` (`run_stage2`) and `src/pga-mpi-merge.c`
(`redistribute_seeds_mem`). The last section covers the timing-barrier fix.

## 1. The problem

The pipeline moves data between ranks twice:

| Exchange | What moves | Record size | Caller |
|---|---|---|---|
| Distributed sort (stage 2) | k-mer records, routed by bucket owner | `rsize` (18 bytes on human chr6) | `run_stage2` in `dsort.c` |
| Seed redistribution | seeds, routed by owner of the contig pair | `1 + IBYTE + JBYTE` | `redistribute_seeds_mem` in `pga-mpi-merge.c` |

Both are *personalized all-to-all*: rank `r` has a different block of bytes
for every peer `p`. MPI's call for that is

```c
MPI_Alltoallv(sendbuf, int sendcounts[], int sdispls[], MPI_BYTE,
              recvbuf, int recvcounts[], int rdispls[], MPI_BYTE, comm);
```

Counts and displacements are `int`, so the maximum is 2^31 − 1 bytes, about
2 GiB. The displacement of the last peer equals the total bytes before it,
so the limit applies to a rank's **total** send and receive volume, not just
to one peer. Above it the values overflow and MPI aborts (or worse,
silently misbehaves).

This hit even at one rank, because stage 2 always exchanges (to itself).
primates4.chr6 has ~8 GB of records and ~23 GB of seeds, so it could not run
at small rank counts. Human chr6 has 11.4 billion records × 18 bytes ≈
206 GB in total, about 3.2 GB per rank at 64 ranks, so it also needs the fix.

Options considered:

- **MPI-4 `MPI_Alltoallv_c`** (64-bit `MPI_Count` arguments). It works on
  Cray MPICH 9.1 (we probed it), but needs MPI-4, and one giant call gives
  no control over staging memory.
- **Bounded rounds of the MPI-3 call** (chosen). Works anywhere, staging is
  bounded, and the result can be made byte-identical.

## 2. The API

```c
void xchg_exchange_sizes(const int64_t *send_bytes, int64_t *recv_bytes, MPI_Comm comm);

void chunked_alltoallv(const int64_t *send_bytes, const int64_t *recv_bytes,
                       int unit, int64_t budget,
                       const XchgSource *src, const XchgSink *dst,
                       MPI_Comm comm, XchgStats *stats);
```

- `send_bytes[p]`: bytes I send to peer `p` (int64).
- `recv_bytes[p]`: bytes `p` sends to me. Obtained with
  `xchg_exchange_sizes`, which is one `MPI_Alltoall` of int64 values (each
  rank sends one number to each peer; that is small and safe).
- `unit`: the record size. Every size must be a multiple of it, so records are
  never split.
- `budget`: per-round bytes, from `xchg_round_budget_from_env()`:
  `PGA_EXCHANGE_ROUND_BYTES` if set and positive, else 1 GiB.
- `src` / `dst`: where data comes from and goes to. Each is either a
  contiguous buffer, or a callback:

```c
typedef struct { const uint8_t *buf; xchg_fill_fn fill; void *ctx; } XchgSource;
typedef struct { uint8_t *buf; xchg_drain_fn drain; void *ctx; } XchgSink;
```

  With `buf`, peer `p`'s data starts at `prefix[p]` (sum of sizes of peers
  before `p`), the standard Alltoallv layout. With a callback,
  `fill(ctx, peer, offset, len, dst)` must copy `len` bytes of that peer's
  data starting at `offset` into `dst`. This lets a caller keep its data in
  separate per-peer buffers instead of packing it into one big array.
- `stats`: rounds used, whether the direct path ran, cap, staging bytes,
  totals. These are the `stage2_rounds`, `redist_n_rounds`, ... fields in the
  `[pr1]` log lines.

## 3. Step by step through `chunked_alltoallv`

### 3.1 Agree on one budget

```c
MPI_Allreduce(MPI_IN_PLACE, &budget, 1, MPI_INT64_T, MPI_MIN, comm);
```

Every rank must make the same decisions (direct vs rounds, cap, number of
rounds), or one rank would call `MPI_Alltoallv` a different number of times
than another and the job would hang. If the environment variable differed
between ranks, taking the minimum makes them agree.

### 3.2 Prefix sums and checks

For each peer it computes `sprefix[p]`, `rprefix[p]` (int64 start offsets in
the full send / receive layout) and the totals `stot`, `rtot`. Any negative
size or size not divisible by `unit` aborts with a message: that would mean a
bug in the caller.

### 3.3 Direct path (the old behaviour)

```c
mine = max(stot, rtot);  Allreduce(MAX) -> worst
contiguous = (src->buf && dst->buf);  Allreduce(MIN) -> all_contig
if (all_contig && worst <= budget && worst <= INT_MAX)  -> one MPI_Alltoallv
```

If **every** rank's totals fit in the budget and everyone uses plain buffers,
it does exactly the old single call. Small datasets therefore run exactly as
before PR #10. The check is global (Allreduce) for the same reason as 3.1:
all ranks must choose the same path.

Note that seed redistribution passes a callback as its source, so
`all_contig` is 0 and it **always** takes the round path, often with
`rounds = 1`. That is why the human logs show `redist_n_rounds=1` with
non-zero staging.

### 3.4 Chunk size `cap`

```c
int64_t xchg_chunk_cap(int npeers, int unit, int64_t budget)
{
    cap  = (budget / npeers) / unit * unit;          // share of budget, whole records
    hard = ((int64_t) INT_MAX / npeers) / unit * unit; // int-safe ceiling
    if (cap > hard) cap = hard;
    if (cap < unit) cap = unit;                       // always make progress
    return cap;
}
```

`cap` is the most a rank sends to **one** peer in **one** round. Why these
bounds:

- `budget / npeers`: a round sends at most `cap` to each of `npeers` peers,
  so the round's send staging is at most `npeers × cap ≤ budget`. Same for
  receive.
- `/ unit * unit`: round down to whole records, so a record never straddles
  two rounds.
- `INT_MAX / npeers`: within one round, the packed displacements reach up to
  `npeers × cap`, which must fit in an `int`.
- `≥ unit`: with a tiny budget and many peers, `cap` would round to 0 and the
  loop would never finish; one record per peer per round always progresses.

Example, human chr6 at 64 ranks: budget 2^30, `rsize` 18.
`2^30 / 64 = 16,777,216`; rounded down to a multiple of 18 gives
`cap = 16,777,206` bytes (932,067 records) per peer per round. The INT_MAX
bound is `33,554,431` → `33,554,430`, larger, so it doesn't bind. Check:
`64 × 16,777,206 × 2 = 2,147,482,368`, exactly the logged `stage2_staging`.

### 3.5 Number of rounds

```c
local_rounds = max over p of ceil(send_bytes[p] / cap), ceil(recv_bytes[p] / cap)
Allreduce(MAX) -> rounds
```

All ranks run `rounds` iterations, the worst case over the whole job. A peer
whose data is already done contributes 0 bytes in later rounds; MPI is fine
with zero counts. Running the same count everywhere keeps the collective
calls matched.

On human chr6 that gave 5 rounds at 64 ranks and 8 at 32 (logged as
`stage2_rounds`). More ranks means a smaller `cap` per peer but also less
data per peer, so the count doesn't simply scale.

### 3.6 Staging buffers

```c
sstage = sum over p of min(send_bytes[p], cap)
rstage = sum over p of min(recv_bytes[p], cap)
```

These are the largest amount any single round can need, so one allocation of
each is reused for every round. Both are ≤ budget, so staging ≤ 2 × budget.
The logged `stage2_staging=2147482368` on human (64 ranks) is just under
2 GiB, which is that bound.

Important: staging is the only memory the budget controls. The caller still
holds the full send data and the full receive buffer (see report §3.2 and
§4.3). The budget limits the *extra* copies, not peak memory.

### 3.7 One round

For round `k` (0-based):

```c
s = xchg_plan_round(send_bytes, n, unit, cap, k, scounts, sdispls);
r = xchg_plan_round(recv_bytes, n, unit, cap, k, rcounts, rdispls);
```

`xchg_plan_round` fills plain `int` arrays for this round only:

- `counts[p] = xchg_round_count(bytes[p], cap, k)`, which is `cap` while the
  peer still has at least `cap` bytes left, the remainder in its last round,
  and 0 afterwards.
- `displs[p]` = running sum of the counts, i.e. packed offsets in the
  staging buffer.
- It returns −1/−2/−3 for negative size, size not a multiple of `unit`, or a
  round total above INT_MAX (which `cap` already prevents).

Then:

1. **Pack**: for every peer with data this round, copy its chunk into
   `sbuf + sdispls[p]`. The chunk starts at offset `k × cap` in that peer's
   range: from `src->buf + sprefix[p] + k×cap`, or via `fill(..., p, k×cap, ...)`.
2. **Exchange**: one ordinary `MPI_Alltoallv` on the staging buffers. All
   counts and displacements fit in `int` by construction.
3. **Unpack**: for every peer with data, copy from `rbuf + rdispls[p]` to
   `dst->buf + rprefix[p] + k×cap`, or call `drain`.

### 3.8 Why the result is byte-identical

In a single Alltoallv, the bytes from peer `p` land at
`[rprefix[p], rprefix[p] + recv_bytes[p])` in their original order. With
rounds, round `k` delivers bytes `[k×cap, k×cap + count)` of that same range
to exactly `rprefix[p] + k×cap`. The rounds cover the range once, in order,
with no gaps or overlaps. So the final receive buffer is the same byte for
byte, for any `cap` and any number of rounds. That is why the
end-to-end tests can demand identical per-rank PAF files rather than just
matching alignment counts.

## 4. The two callers

### 4.1 Stage 2 (`run_stage2`, `dsort.c`)

1. `organize_send_buffer` packs the local records into one buffer grouped by
   destination rank, giving `send_counts` (bytes per peer, int64).
2. `xchg_exchange_sizes` → `recv_counts`.
3. Allocate `recv->data` with `recv_total + 1` bytes (the +1 is the
   `msd_sort` sentinel; keep it).
4. `chunked_alltoallv(..., rsize, budget, {send_buf}, {recv->data}, ...)`.
   Both sides are plain buffers, so small runs take the direct path.
5. Free the send buffer, then sort locally and fix LCPs.

During the call a rank holds its extracted records (`local`), the packed copy
(`send_buf`), the received slice, and the staging buffers. That's why stage 2
is often the memory peak.

### 4.2 Seed redistribution (`redistribute_seeds_mem`, `pga-mpi-merge.c`)

The merge writes each seed straight into a growing buffer for its
destination rank (`SeedUnit`, one per destination, separately for same-strand
"N" and opposite-strand "C" seeds). Then:

1. `send_sizes[r] = units[r].len`; exchange sizes.
2. Allocate one receive buffer of the total.
3. Source = callback `seed_mem_fill`, which just does
   `memcpy(dst, units[peer].data + off, len)`. No packing into one big send
   array: each round copies straight from the per-destination buffers into
   staging.
4. After the call, free all send units.
5. On one rank the function skips MPI entirely and hands over the buffer.

The received buffer is in source-rank order, and each source's seeds are in
its merge order. That is global k-mer order, independent of the rank count.

## 5. Tests and validation

| Level | What | Where |
|---|---|---|
| Planner, no MPI | `xchg_chunk_cap`, `xchg_rounds_for`, `xchg_round_count`, `xchg_plan_round` with sizes far above 2^31, no allocation | `tests/unit/test_chunked_exchange.c`, `test_exchange_plan.c` (R-02, R-06) |
| MPI | Exchanges vs a reference on 1/2/3/4/8 ranks, default budget and `PGA_EXCHANGE_ROUND_BYTES=256` (many rounds) | `tests/mpi/test_chunked_mpi.c`, `make test-mpi` |
| Phase verifier | yeast7 `verify-phase3-stage2`, 1 MiB budget, N = 1/2/4/8 | `docs/agent/pr2-results.md` |
| End to end | yeast7 N=1/2/8, scerevisiae8 N=4 at 8 MiB budget: raw per-rank PAF byte-identical to the pre-O2 build, gold MATCH | `docs/agent/pr2-results.md` |

Wall time was unchanged within noise even with 100+ rounds (yeast7 N=1: 64.8
vs 65.5 s).

## 6. Things to know when changing it

- **Every decision must be collective.** Budget, path, cap and round count
  are all agreed with Allreduce. Any per-rank shortcut (e.g. "I have nothing
  to send, skip the loop") deadlocks the others.
- **Sizes must be multiples of `unit`.** Callers guarantee it; the routine
  aborts otherwise.
- **The budget is not a memory cap** for the whole exchange, only for
  staging. Lowering it adds rounds without reducing the full send/receive
  buffers.
- **The seed path always goes through rounds** because its source is a
  callback. Making it use the direct path would need a contiguous send
  buffer, which is the copy O4 removed.
- `XchgStats` feed the `[pr1]` log line; scripts read
  `stage2_rounds`, `redist_n_rounds`, `redist_c_rounds`, `*_staging`.

## 7. Timing-barrier fix (branch `timing-barrier`)

Raised by Saksham (2026-11-14). `pga-mpi-merge.c` timed the extract stage
without synchronizing. A rank that finished extraction early went straight
into `run_stage2`, whose first call is the histogram `MPI_Allreduce`, and
waited there for the slowest extractor. The stage table reports the max
over ranks of each stage, so:

- `extract` = slowest rank's extraction (correct), and
- `sort` = some early rank's sort, **including** its wait for that slowest
  extraction.

The wait was counted twice, and `TOTAL` (sum of the four maxima) came out
above the real wall time. It also inflated the per-rank `hist_s` sub-timer:
on human chr6 most of its ~2.4–3 s is that wait.

Fix:

```c
extract_to_records(gdb, world_size, my_rank, &local);
double t_extract_own = MPI_Wtime() - t_stage_start;   // this rank's own work
MPI_Barrier(MPI_COMM_WORLD);
t_extract = MPI_Wtime() - t_stage_start;              // stage time, all ranks done
```

- The stage table and CSV use `t_extract` (after the barrier), so the sort
  starts synchronized and nothing is counted twice.
- The per-rank `[o1-0]` line keeps `extract_s` = the rank's own time (same
  meaning as before, so O1 imbalance analyses stay comparable) and adds
  `extract_wait_s`.

The other stage boundaries are already synchronized. The sort ends with
the cross-rank boundary check (`MPI_Allgather` + `MPI_Allreduce`), and the
merge stage ends with the seed exchange, so no further barriers are
needed. The barrier changes no output, only timings. Timings from runs
before this fix overstate `sort_s` and `TOTAL` by up to the extraction
imbalance.
