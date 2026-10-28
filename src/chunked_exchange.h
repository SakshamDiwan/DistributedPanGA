#ifndef PGA_CHUNKED_EXCHANGE_H
#define PGA_CHUNKED_EXCHANGE_H

// All-to-all byte exchange whose sizes may exceed INT_MAX, done as bounded
// rounds of MPI-3 MPI_Alltoallv.
//
// Every peer's data is a contiguous byte range of a known int64 size. Each
// round moves at most `cap` bytes per peer, and `cap` is a multiple of the
// record size, so no record is split across rounds. Received chunks are
// placed at their final offsets, so the result is byte-identical to one
// unbounded Alltoallv.
//
// The planner functions below are MPI-free so they can be tested with huge
// synthetic sizes and no allocation.

#include <stdint.h>
#include <mpi.h>

// Per-round budget: PGA_EXCHANGE_ROUND_BYTES if set (> 0), else this.
#define XCHG_DEFAULT_ROUND_BYTES ((int64_t) 1 << 30)
int64_t xchg_round_budget_from_env(void);

// Bytes per peer per round: floor(budget / npeers) rounded down to a multiple
// of unit, at least unit, and at most INT_MAX / npeers (rounded down to unit),
// so a round's packed displacements always fit in an int.
int64_t xchg_chunk_cap(int npeers, int unit, int64_t budget);

// Rounds needed to move `bytes` in chunks of `cap`: ceil(bytes / cap).
int64_t xchg_rounds_for(int64_t bytes, int64_t cap);

// Bytes moved for one peer in round k (0 once that peer's data is done).
int64_t xchg_round_count(int64_t bytes, int64_t cap, int64_t k);

// Fill the int counts[]/displs[] for round k (displs are packed prefix sums
// into the round's staging buffer). Returns the round's total bytes, or
// -1 if a size is negative, -2 if a size is not a multiple of unit,
// -3 if a round total would exceed INT_MAX (cannot happen with xchg_chunk_cap).
int64_t xchg_plan_round(const int64_t *bytes, int npeers, int unit, int64_t cap,
                        int64_t k, int *counts, int *displs);

// Data access. Either buf != NULL (peer p's bytes start at prefix[p] =
// sum of sizes of peers < p) or the callback is used. offset is relative to
// the start of that peer's range.
typedef void (*xchg_fill_fn)(void *ctx, int peer, int64_t offset, int64_t len,
                             uint8_t *dst);
typedef void (*xchg_drain_fn)(void *ctx, int peer, int64_t offset, int64_t len,
                              const uint8_t *src);
typedef struct { const uint8_t *buf; xchg_fill_fn fill; void *ctx; } XchgSource;
typedef struct { uint8_t *buf; xchg_drain_fn drain; void *ctx; } XchgSink;

typedef struct {
    int64_t rounds;          // 1 on the direct path
    int     direct;          // 1 if a single MPI_Alltoallv was used
    int64_t cap;             // bytes per peer per round (chunked path)
    int64_t staging_bytes;   // send + receive staging allocated (chunked path)
    int64_t sent_bytes;
    int64_t recv_bytes;
} XchgStats;

// Exchange the int64 per-peer sizes: recv_bytes[p] = what p sends to me.
void xchg_exchange_sizes(const int64_t *send_bytes, int64_t *recv_bytes, MPI_Comm comm);

// Collective. send_bytes/recv_bytes as from xchg_exchange_sizes; every size
// must be a multiple of unit. The ranks use the minimum of their budgets. Uses one MPI_Alltoallv when every rank's send
// and receive totals fit in the budget and both sides are contiguous buffers;
// otherwise rounds. Aborts the job on allocation or planning failure.
void chunked_alltoallv(const int64_t *send_bytes, const int64_t *recv_bytes,
                       int unit, int64_t budget,
                       const XchgSource *src, const XchgSink *dst,
                       MPI_Comm comm, XchgStats *stats);

#endif
