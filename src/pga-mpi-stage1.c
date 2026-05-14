// pga-mpi-stage1: Phase 3 stage 1 — gather all records to rank 0, then sort.
// Throwaway binary; gets deleted once stage 2 is green.
//
// Pipeline:
//   1. Every rank: Read_GDB, compute_work_split(world_size, my_rank), extract_to_records.
//   2. mpi_gather_records to rank 0.
//   3. Rank 0: sort_records, in-buffer order check, dump.
//
// Usage:  mpirun -n N pga-mpi-stage1 <gdb-stem> [output-tuples-file]

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>

#include "extract.h"
#include "record.h"
#include "sort.h"
#include "dsort.h"
#include "constants.h"
#include "GDB.h"

#define ROOT 0

static const char hex_chars[] = "0123456789abcdef";

static int64_t read_le_field(const uint8_t *bytes, int n)
{
    int64_t v = 0;
    for (int k = 0; k < n; k++)
        v |= ((int64_t) bytes[k]) << (8 * k);
    return v;
}

static void dump_records_as_tuples(FILE *out, const RecordBuffer *buf)
{
    int rsize       = buf->sizing.record_size;
    int post_bytes  = buf->sizing.post_bytes;
    int cont_bytes  = buf->sizing.cont_bytes;
    int pos_off     = POS_OFFSET;
    int cont_off    = POS_OFFSET + post_bytes;

    int64_t strand_mask = ((int64_t) 1) << (8 * cont_bytes - 1);
    int64_t contig_mask = strand_mask - 1;

    char line[128];
    for (int64_t i = 0; i < buf->count; i++)
    {
        const uint8_t *r = buf->data + i * rsize;

        int64_t pos     = read_le_field(r + pos_off, post_bytes);
        int64_t cont_v  = read_le_field(r + cont_off, cont_bytes);
        int     strand  = (cont_v & strand_mask) ? 1 : 0;
        int     contig  = (int) (cont_v & contig_mask);

        int n = snprintf(line, sizeof(line), "%d %lld %d ",
                         contig, (long long) pos, strand);
        for (int b = 0; b < KBYTES; b++)
        {
            uint8_t kb = r[KMER_OFFSET + b];
            line[n++] = hex_chars[kb >> 4];
            line[n++] = hex_chars[kb & 0xf];
        }
        line[n++] = '\n';
        fwrite(line, 1, n, out);
    }
}

// Same in-buffer order check we added to pga-sort: walk records in place,
// assert non-decreasing k-mer bytes. Returns 0 on PASS, nonzero violation count on FAIL.
static int64_t check_sort_order(const RecordBuffer *buf, int64_t *first_bad_out)
{
    int rsize = buf->sizing.record_size;
    int64_t bad = 0;
    int64_t first_bad = -1;
    for (int64_t i = 1; i < buf->count; i++)
    {
        const uint8_t *a = buf->data + (i - 1) * rsize + KMER_OFFSET;
        const uint8_t *b = buf->data +  i      * rsize + KMER_OFFSET;
        if (memcmp(a, b, KBYTES) > 0)
        {
            if (first_bad < 0) first_bad = i;
            bad++;
        }
    }
    if (first_bad_out) *first_bad_out = first_bad;
    return bad;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc < 2 || argc > 3)
    {
        if (my_rank == ROOT)
            fprintf(stderr, "Usage: %s <gdb-stem> [output-tuples-file]\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    GDB _gdb, *gdb = &_gdb;
    if (Read_GDB(gdb, argv[1]) < 0)
    {
        fprintf(stderr, "rank %d: Read_GDB failed for '%s'\n", my_rank, argv[1]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (my_rank == ROOT)
        fprintf(stderr, "[stage1] world_size=%d  GDB: %d contigs, %lld total bases\n",
                world_size, gdb->ncontig, (long long) gdb->seqtot);

    RecordSizing sizing = compute_record_sizing(gdb);
    if (my_rank == ROOT)
        fprintf(stderr, "[stage1] record layout: post=%d cont=%d size=%d\n",
                sizing.post_bytes, sizing.cont_bytes, sizing.record_size);

    RecordBuffer local;
    record_buffer_init(&local, sizing);
    extract_to_records(gdb, world_size, my_rank, &local);

    fprintf(stderr, "[stage1] rank %d: extracted %lld records\n",
            my_rank, (long long) local.count);

    RecordBuffer gathered;
    memset(&gathered, 0, sizeof(gathered));
    mpi_gather_records(&local, &gathered, ROOT, MPI_COMM_WORLD);

    record_buffer_free(&local);

    if (my_rank == ROOT)
    {
        fprintf(stderr, "[stage1] gathered %lld records on root, sorting...\n",
                (long long) gathered.count);

        sort_records(&gathered);

        int64_t first_bad = -1;
        int64_t bad = check_sort_order(&gathered, &first_bad);
        if (bad == 0)
            fprintf(stderr, "[stage1] sort order check: PASS (%lld records)\n",
                    (long long) gathered.count);
        else
        {
            fprintf(stderr, "[stage1] sort order check: FAIL (%lld violations, first at %lld)\n",
                    (long long) bad, (long long) first_bad);
            record_buffer_free(&gathered);
            Close_GDB(gdb);
            MPI_Abort(MPI_COMM_WORLD, 2);
        }

        FILE *out = stdout;
        if (argc == 3)
        {
            out = fopen(argv[2], "w");
            if (out == NULL)
            {
                fprintf(stderr, "[stage1] cannot open '%s' for write\n", argv[2]);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        fprintf(stderr, "[stage1] writing tuples...\n");
        dump_records_as_tuples(out, &gathered);
        if (out != stdout) fclose(out);

        record_buffer_free(&gathered);
    }

    Close_GDB(gdb);
    MPI_Finalize();
    return 0;
}
