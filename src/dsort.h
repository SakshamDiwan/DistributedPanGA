#ifndef DPGA_DSORT_H
#define DPGA_DSORT_H

#include <stdint.h>
#include <mpi.h>
#include "extract.h"
#include "constants.h"

// Stage 1 — gather to rank 0.
// Every rank's local RecordBuffer is gathered onto `root`. On the root, `out`
// is initialized to hold all records (with the same RecordSizing as `local`);
// on non-roots, `out` is left untouched.
//
// Both buffers must use the same RecordSizing (same record_size, post_bytes,
// cont_bytes). Records are sent as MPI_BYTE; no custom datatype.
//
// On return:
//   - root rank: out->count = sum of all ranks' local->count, out->data
//                holds rank-0 records first, then rank-1, etc.
//   - non-root:  out is untouched (caller should not have init'd it).
//
// The `local` buffer is NOT freed by this call.
void mpi_gather_records(const RecordBuffer *local, RecordBuffer *out,
                        int root, MPI_Comm comm);


// Stage 2 — Alltoallv redistribution + local sort + LCP fixup.

// bucket = first 5 bases = top byte | top 2 bits of next byte.
static inline int record_bucket(const uint8_t *kmer)
{
    return ((int) kmer[0] << 2) | ((int) kmer[1] >> 6);
}

// Walk `local` and count records per 10-bit bucket. Caller-allocated
// `local_counts[NUM_BUCK]` is overwritten.
void bucket_local_counts(const RecordBuffer *local, int64_t *local_counts);

// Compute the destination map from a global histogram.
//   global_counts[NUM_BUCK]  — input, summed across all ranks
//   ksplit[world_size+1]     — output, bucket-index split points
//   select[NUM_BUCK]         — output, Select[bucket] = destination rank
// Algorithm mirrors FastGA's distribute() lines 678-697 with NPARTS = world_size.
// Deterministic given identical input → safe to compute redundantly on every rank.
void compute_destination_map(const int64_t *global_counts, int world_size,
                             int *ksplit, int *select);

// Reorganize `local` into a contiguous send buffer ordered by destination rank.
// Out-of-place per-rank lists, then concat.
//
//   local            — input buffer (read-only; caller frees later)
//   select           — Select[bucket] -> dest rank, length NUM_BUCK
//   world_size       — number of destination slots
//   send_buf_out     — *send_buf_out set to a freshly malloc'd buffer of
//                      total local->count records, grouped by destination
//   send_counts_out  — caller-supplied int[world_size], filled with byte counts
//   send_displs_out  — caller-supplied int[world_size], filled with byte displs
//
// Caller frees *send_buf_out when done.
void organize_send_buffer(const RecordBuffer *local, const int *select,
                          int world_size,
                          uint8_t **send_buf_out,
                          int *send_counts_out, int *send_displs_out);

// Run the full stage 2 pipeline on this rank.
// Inputs:
//   local      — this rank's extracted records (from extract_to_records)
//   comm       — MPI communicator (typically MPI_COMM_WORLD)
// Output:
//   recv       — initialized to hold this rank's slice of redistributed records,
//                sorted by k-mer with LCP byte filled correctly relative to the
//                global predecessor (cross-rank LCP boundary fixed up).
//
// Caller must free `recv` (record_buffer_free) when done.
// `local` is NOT freed by this call.
void run_stage2(const RecordBuffer *local, RecordBuffer *recv, MPI_Comm comm);


// Walk recv->data, assert non-decreasing k-mer-byte order.
// Returns: number of out-of-order pairs (0 = sorted). If first_bad is non-NULL
// and a violation is found, *first_bad is set to the offending record index.
int64_t verify_per_rank_sorted(const RecordBuffer *recv, int64_t *first_bad);

// Cross-rank check: every rank r exchanges its boundary k-mers with rank r-1
// and rank r+1. Asserts last_kmer[r] <= first_kmer[r+1]. Also returns this
// rank's predecessor's last k-mer in `prev_last_kmer_out` (size KBYTES) so the
// LCP fixup can reuse it. For rank 0, prev_last_kmer_out is left zeroed.
//
// Returns: 0 if all boundaries are ordered, nonzero if this rank's first k-mer
// is < the previous rank's last k-mer.
// `have_prev_out` (may be NULL) reports whether a NON-EMPTY predecessor rank
// exists at all. It is 0 for rank 0, and also for any rank whose predecessors
// all hold zero records -- in which case this rank owns the globally first
// record and its LCP must stay 0 rather than being measured against anything.
int verify_cross_rank_boundaries(const RecordBuffer *recv,
                                 uint8_t *prev_last_kmer_out,
                                 int *have_prev_out,
                                 MPI_Comm comm);

#endif
