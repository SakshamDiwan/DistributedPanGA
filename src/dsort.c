#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <mpi.h>

#include "dsort.h"
#include "sort.h"
#include "record.h"

// LCP encoding.
//
// Byte 0 of every record holds a SEMANTIC LCP — the number of bases in common
// with the previous record's k-mer, range [0, 40]. This is computed by
// simple_kmer_lcp_bases (in sort.c) and applied uniformly across the buffer by
// recalc_all_lcps after sort_records, overwriting whatever msd_sort wrote.
//
// We deliberately do NOT keep msd_sort's internal byte 0 because its encoding
// mixes two formulas (intra-bucket and inter-bucket) that don't agree with
// "bases in common" at part[] boundaries. Same convention applies at rank
// seams (the cross-rank fixup in run_stage2 also writes simple_kmer_lcp_bases).

// Stage 1 — gather to rank 0.

void mpi_gather_records(const RecordBuffer *local, RecordBuffer *out,
                        int root, MPI_Comm comm)
{
    int my_rank, world_size;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &world_size);

    int rsize = local->sizing.record_size;

    int64_t local_bytes_64 = local->count * (int64_t) rsize;
    if (local_bytes_64 > INT_MAX)
    {
        fprintf(stderr, "rank %d: local byte count %lld exceeds int range\n",
                my_rank, (long long) local_bytes_64);
        MPI_Abort(comm, 1);
    }
    int local_bytes = (int) local_bytes_64;

    int *recv_counts = NULL;
    int *recv_displs = NULL;
    if (my_rank == root)
    {
        recv_counts = (int *) malloc(world_size * sizeof(int));
        recv_displs = (int *) malloc(world_size * sizeof(int));
        if (recv_counts == NULL || recv_displs == NULL)
        {
            fprintf(stderr, "mpi_gather_records: out of memory for count arrays\n");
            MPI_Abort(comm, 1);
        }
    }

    MPI_Gather(&local_bytes, 1, MPI_INT,
               recv_counts,  1, MPI_INT,
               root, comm);

    if (my_rank == root)
    {
        int64_t total_bytes = 0;
        for (int r = 0; r < world_size; r++)
        {
            recv_displs[r] = (int) total_bytes;
            total_bytes += recv_counts[r];
            if (total_bytes > INT_MAX)
            {
                fprintf(stderr,
                        "mpi_gather_records: total byte count %lld exceeds int range "
                        "(MPI_Gatherv displ uses int)\n",
                        (long long) total_bytes);
                MPI_Abort(comm, 1);
            }
        }
        int64_t total_records = total_bytes / rsize;

        out->sizing   = local->sizing;
        out->count    = total_records;
        out->capacity = total_records;
        out->data = (uint8_t *) malloc(total_bytes + 1);
        if (out->data == NULL)
        {
            fprintf(stderr,
                    "mpi_gather_records: out of memory for gathered buffer (%lld bytes)\n",
                    (long long) total_bytes);
            MPI_Abort(comm, 1);
        }
    }

    MPI_Gatherv(local->data, local_bytes, MPI_BYTE,
                (my_rank == root) ? out->data : NULL,
                recv_counts, recv_displs, MPI_BYTE,
                root, comm);

    if (my_rank == root)
    {
        free(recv_counts);
        free(recv_displs);
    }
}


// Stage 2 — bucket counting and splitter.

void bucket_local_counts(const RecordBuffer *local, int64_t *local_counts)
{
    int rsize = local->sizing.record_size;
    for (int i = 0; i < NUM_BUCK; i++)
        local_counts[i] = 0;
    for (int64_t i = 0; i < local->count; i++)
    {
        const uint8_t *kmer = local->data + i * rsize + KMER_OFFSET;
        local_counts[record_bucket(kmer)] += 1;
    }
}

void compute_destination_map(const int64_t *global_counts, int world_size,
                             int *ksplit, int *select)
{
    // Mirror of FastGA's distribute() lines 678-697 with NPARTS = world_size.
    int64_t cum[NUM_BUCK];
    cum[0] = global_counts[0];
    for (int i = 1; i < NUM_BUCK; i++)
        cum[i] = cum[i-1] + global_counts[i];

    int64_t total = cum[NUM_BUCK - 1];

    ksplit[0] = 0;
    int n = 1;
    int64_t t = total / world_size;

    for (int i = 0; i < NUM_BUCK; i++)
    {
        if (cum[i] >= t)
        {
            int64_t prev = (i > 0) ? cum[i-1] : 0;
            if (cum[i] - t > t - prev)
            {
                select[i]   = n;
                ksplit[n]   = i;
            }
            else
            {
                select[i]   = n - 1;
                ksplit[n]   = i + 1;
            }
            n += 1;
            // Guard: don't overshoot the array if rounding produces too many splits.
            if (n > world_size) n = world_size;
            t = ((int64_t) n * total) / world_size;
        }
        else
        {
            select[i] = n - 1;
        }
    }
    ksplit[world_size] = NUM_BUCK;

    // Defensive: if rounding/empty buckets produced fewer splits than world_size,
    // fill the trailing ksplit entries with NUM_BUCK so the empty ranks get
    // empty (but well-defined) slices.
    for (int r = n; r < world_size; r++)
        ksplit[r] = NUM_BUCK;
}


// Stage 2 — send-side organization (per-rank lists, then concat).

void organize_send_buffer(const RecordBuffer *local, const int *select,
                          int world_size,
                          uint8_t **send_buf_out,
                          int *send_counts_out, int *send_displs_out)
{
    int rsize = local->sizing.record_size;

    int64_t per_rank_records[world_size];
    for (int r = 0; r < world_size; r++)
        per_rank_records[r] = 0;

    for (int64_t i = 0; i < local->count; i++)
    {
        const uint8_t *kmer = local->data + i * rsize + KMER_OFFSET;
        int dest = select[record_bucket(kmer)];
        per_rank_records[dest] += 1;
    }

    int64_t total_bytes = 0;
    for (int r = 0; r < world_size; r++)
    {
        int64_t bytes_for_r = per_rank_records[r] * (int64_t) rsize;
        if (bytes_for_r > INT_MAX || total_bytes + bytes_for_r > INT_MAX)
        {
            fprintf(stderr,
                    "organize_send_buffer: per-rank or total send count exceeds int range\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        send_counts_out[r] = (int) bytes_for_r;
        send_displs_out[r] = (int) total_bytes;
        total_bytes += bytes_for_r;
    }

    // Allocate the contiguous send buffer (+1 for any downstream sentinel).
    uint8_t *send_buf = (uint8_t *) malloc(total_bytes + 1);
    if (send_buf == NULL)
    {
        fprintf(stderr,
                "organize_send_buffer: out of memory for send buffer (%lld bytes)\n",
                (long long) total_bytes);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Pass 2: walk records again, write each to its destination's slot.
    // cursor[r] tracks byte offset within send_buf for rank r.
    int64_t cursor[world_size];
    for (int r = 0; r < world_size; r++)
        cursor[r] = (int64_t) send_displs_out[r];

    for (int64_t i = 0; i < local->count; i++)
    {
        const uint8_t *src = local->data + i * rsize;
        int dest = select[record_bucket(src + KMER_OFFSET)];
        memcpy(send_buf + cursor[dest], src, rsize);
        cursor[dest] += rsize;
    }

    *send_buf_out = send_buf;
}


// Stage 2 — verification helpers.

int64_t verify_per_rank_sorted(const RecordBuffer *recv, int64_t *first_bad)
{
    int rsize = recv->sizing.record_size;
    int64_t bad = 0;
    int64_t fb = -1;
    for (int64_t i = 1; i < recv->count; i++)
    {
        const uint8_t *a = recv->data + (i - 1) * rsize + KMER_OFFSET;
        const uint8_t *b = recv->data +  i      * rsize + KMER_OFFSET;
        if (memcmp(a, b, KBYTES) > 0)
        {
            if (fb < 0) fb = i;
            bad++;
        }
    }
    if (first_bad) *first_bad = fb;
    return bad;
}

int verify_cross_rank_boundaries(const RecordBuffer *recv,
                                 uint8_t *prev_last_kmer_out,
                                 int *have_prev_out,
                                 MPI_Comm comm)
{
    int have_prev = 0;
    int my_rank, world_size;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &world_size);

    int rsize = recv->sizing.record_size;

    // Edge k-mers we own. If we have zero records, send a max-valued sentinel
    // up to our successor (so its boundary check passes vacuously) and ignore
    // what we receive from our predecessor for the comparison.
    uint8_t my_first[KBYTES];
    uint8_t my_last[KBYTES];
    int     have_records = (recv->count > 0);

    if (have_records)
    {
        memcpy(my_first, recv->data + KMER_OFFSET, KBYTES);
        memcpy(my_last,
               recv->data + (recv->count - 1) * rsize + KMER_OFFSET, KBYTES);
    }
    else
    {
        memset(my_first, 0xff, KBYTES);  // empty: any predecessor will be <= this
        memset(my_last,  0x00, KBYTES);  // empty: any successor will be >= this
    }

    uint8_t prev_last[KBYTES];
    memset(prev_last, 0, KBYTES);

    // Find the last k-mer of the nearest NON-EMPTY rank before this one.
    //
    // A one-hop MPI_Sendrecv is wrong here: a rank holding no records has no
    // last k-mer, so it would contribute an all-zero sentinel and its successor
    // would compute its seam LCP against that instead of against the real
    // predecessor, which may be several ranks back. The boundary check would
    // also pass vacuously against the sentinel.
    //
    // World size is bounded (contig assignment caps it at 127) and the payload
    // is KBYTES+1 per rank, so an Allgather is cheap and obviously correct.
    {
        int      stride = KBYTES + 1;
        uint8_t  mine[KBYTES + 1];
        uint8_t *all = (uint8_t *) malloc((size_t) world_size * stride);

        if (all == NULL)
        {
            fprintf(stderr, "verify_cross_rank_boundaries: out of memory\n");
            MPI_Abort(comm, 1);
        }

        mine[0] = (uint8_t) (have_records ? 1 : 0);
        memcpy(mine + 1, my_last, KBYTES);
        MPI_Allgather(mine, stride, MPI_BYTE, all, stride, MPI_BYTE, comm);

        for (int r = my_rank - 1; r >= 0; r--)
            if (all[(size_t) r * stride])                 // that rank had records
            {
                memcpy(prev_last, all + (size_t) r * stride + 1, KBYTES);
                have_prev = 1;
                break;
            }
        free(all);
    }

    int local_bad = 0;
    if (have_prev && have_records)
    {
        if (memcmp(prev_last, my_first, KBYTES) > 0)
            local_bad = 1;
    }

    if (prev_last_kmer_out)
        memcpy(prev_last_kmer_out, prev_last, KBYTES);
    if (have_prev_out)
        *have_prev_out = have_prev;

    int global_bad;
    MPI_Allreduce(&local_bad, &global_bad, 1, MPI_INT, MPI_SUM, comm);
    return global_bad;
}


// Stage 2 — full pipeline.

void run_stage2(const RecordBuffer *local, RecordBuffer *recv, MPI_Comm comm)
{
    int my_rank, world_size;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &world_size);

    int rsize = local->sizing.record_size;

    int64_t local_counts[NUM_BUCK];
    bucket_local_counts(local, local_counts);

    int64_t global_counts[NUM_BUCK];
    MPI_Allreduce(local_counts, global_counts, NUM_BUCK,
                  MPI_LONG_LONG, MPI_SUM, comm);

    // Splitter computed identically on every rank (deterministic given identical input).
    int *ksplit = (int *) malloc((world_size + 1) * sizeof(int));
    int *select = (int *) malloc(NUM_BUCK * sizeof(int));
    if (ksplit == NULL || select == NULL)
    {
        fprintf(stderr, "run_stage2: out of memory for splitter arrays\n");
        MPI_Abort(comm, 1);
    }
    compute_destination_map(global_counts, world_size, ksplit, select);

    int *send_counts = (int *) malloc(world_size * sizeof(int));
    int *send_displs = (int *) malloc(world_size * sizeof(int));
    int *recv_counts = (int *) malloc(world_size * sizeof(int));
    int *recv_displs = (int *) malloc(world_size * sizeof(int));
    if (!send_counts || !send_displs || !recv_counts || !recv_displs)
    {
        fprintf(stderr, "run_stage2: out of memory for count arrays\n");
        MPI_Abort(comm, 1);
    }

    uint8_t *send_buf = NULL;
    organize_send_buffer(local, select, world_size,
                         &send_buf, send_counts, send_displs);

    MPI_Alltoall(send_counts, 1, MPI_INT,
                 recv_counts, 1, MPI_INT, comm);

    int64_t recv_total_bytes = 0;
    for (int r = 0; r < world_size; r++)
    {
        recv_displs[r] = (int) recv_total_bytes;
        recv_total_bytes += recv_counts[r];
        if (recv_total_bytes > INT_MAX)
        {
            fprintf(stderr,
                    "run_stage2: incoming byte count %lld on rank %d exceeds int range\n",
                    (long long) recv_total_bytes, my_rank);
            MPI_Abort(comm, 1);
        }
    }
    int64_t recv_record_count = recv_total_bytes / rsize;

    recv->sizing   = local->sizing;
    recv->count    = recv_record_count;
    recv->capacity = recv_record_count;
    recv->data     = (uint8_t *) malloc(recv_total_bytes + 1);
    if (recv->data == NULL)
    {
        fprintf(stderr, "run_stage2: out of memory for recv buffer on rank %d\n", my_rank);
        MPI_Abort(comm, 1);
    }

    MPI_Alltoallv(send_buf,   send_counts, send_displs, MPI_BYTE,
                  recv->data, recv_counts, recv_displs, MPI_BYTE,
                  comm);

    free(send_buf);
    free(send_counts);
    free(send_displs);
    free(recv_counts);
    free(recv_displs);
    free(ksplit);
    free(select);

    sort_records(recv);

#ifdef LCPs
    // Overwrite msd_sort's mixed-formula byte 0 with our base-level LCP.
    recalc_all_lcps(recv);
#endif
    int64_t first_bad = -1;
    int64_t bad = verify_per_rank_sorted(recv, &first_bad);
    if (bad != 0)
    {
        fprintf(stderr,
                "run_stage2: rank %d local sort check FAIL (%lld violations, first at %lld)\n",
                my_rank, (long long) bad, (long long) first_bad);
        MPI_Abort(comm, 2);
    }

    // Cross-rank boundary check; also collects prev_last for the LCP fixup below.
    uint8_t prev_last_kmer[KBYTES];
    int have_prev = 0;
    int boundary_bad = verify_cross_rank_boundaries(recv, prev_last_kmer,
                                                    &have_prev, comm);
    if (boundary_bad != 0)
    {
        if (my_rank == 0)
            fprintf(stderr,
                    "run_stage2: cross-rank boundary check FAIL (%d ranks reported violation)\n",
                    boundary_bad);
        MPI_Abort(comm, 3);
    }

    // LCP boundary fixup. For ranks ≥ 1, recalc_all_lcps wrote 0 into
    // recv->data[0] because it treats the rank-local first record as having no
    // predecessor. The actual predecessor is on rank r-1, exchanged above.
    // Overwrite recv->data[0] with the correct cross-rank LCP using the same
    // convention as the rest of the buffer:
    //
    //   LCPS=1 build (default):  simple_kmer_lcp_bases(prev_last, my_first)
    //                            range [0, 40], same as recalc_all_lcps used
    //                            for every other record.
    //
    //   LCPS=0 build (testing):  byte 0 is a 0/1 boundary marker.
    // Only when a non-empty predecessor actually exists. If every rank before
    // this one is empty, this rank holds the globally first record, whose LCP
    // recalc_all_lcps already set to 0 -- overwriting it here would measure it
    // against a sentinel that corresponds to no real k-mer.
    if (have_prev && recv->count > 0)
    {
        const uint8_t *my_first = recv->data + KMER_OFFSET;
#ifdef LCPs
        recv->data[0] = (uint8_t) simple_kmer_lcp_bases(prev_last_kmer, my_first, KBYTES);
#else
        int same = (memcmp(prev_last_kmer, my_first, KBYTES) == 0);
        recv->data[0] = same ? 0 : 1;
#endif
    }
}
