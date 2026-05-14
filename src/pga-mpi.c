// pga-mpi: Phase 3 stage 2 — distributed extract + Alltoallv redistribute +
// per-rank sort. Each rank dumps its sorted slice to <out-stem>.<rank>.tuples.
//
// Pipeline:
//   1. Every rank: Read_GDB, compute_work_split(world_size, my_rank), extract.
//   2. run_stage2: bucket count + Allreduce + splitter + Alltoallv + local sort
//      + per-rank order check + cross-rank boundary check + LCP fixup.
//   3. Each rank dumps recv buffer as text tuples.
//
// Usage:  mpirun -n N pga-mpi <gdb-stem> <output-stem>
//   Writes: <output-stem>.0.tuples, <output-stem>.1.tuples, ..., <output-stem>.{N-1}.tuples

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

// Two output formats. The default 4-column format keeps existing verify
// scripts (Tier 1 multiset, Tier 2 order check) unchanged. The 5-column
// format is opt-in via env var PGA_DUMP_LCP=1 and adds byte 0 (the LCP byte
// as msd_sort wrote it, with our cross-rank fixup) as the second-to-last
// column. The k-mer hex stays last so existing awk scripts that key on $4
// for the k-mer continue to work — LCP-aware scripts key on $4 (lcp) and
// $5 (kmer hex).
//
// Default:        "<contig> <position> <strand> <hex_kmer>\n"
// PGA_DUMP_LCP=1: "<contig> <position> <strand> <lcp> <hex_kmer>\n"
static void dump_records_as_tuples(FILE *out, const RecordBuffer *buf, int with_lcp)
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

        int n;
        if (with_lcp)
            n = snprintf(line, sizeof(line), "%d %lld %d %u ",
                         contig, (long long) pos, strand, (unsigned) r[LCP_OFFSET]);
        else
            n = snprintf(line, sizeof(line), "%d %lld %d ",
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

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc != 3)
    {
        if (my_rank == ROOT)
            fprintf(stderr, "Usage: %s <gdb-stem> <output-stem>\n", argv[0]);
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
        fprintf(stderr, "[stage2] world_size=%d  GDB: %d contigs, %lld total bases\n",
                world_size, gdb->ncontig, (long long) gdb->seqtot);

    RecordSizing sizing = compute_record_sizing(gdb);

    RecordBuffer local;
    record_buffer_init(&local, sizing);
    extract_to_records(gdb, world_size, my_rank, &local);

    fprintf(stderr, "[stage2] rank %d: extracted %lld records\n",
            my_rank, (long long) local.count);

    RecordBuffer recv;
    memset(&recv, 0, sizeof(recv));
    run_stage2(&local, &recv, MPI_COMM_WORLD);

    record_buffer_free(&local);

    fprintf(stderr, "[stage2] rank %d: own slice has %lld records (post-redistribute, sorted)\n",
            my_rank, (long long) recv.count);

    char outpath[1024];
    snprintf(outpath, sizeof(outpath), "%s.%d.tuples", argv[2], my_rank);
    FILE *out = fopen(outpath, "w");
    if (out == NULL)
    {
        fprintf(stderr, "rank %d: cannot open '%s' for write\n", my_rank, outpath);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int with_lcp = (getenv("PGA_DUMP_LCP") != NULL && getenv("PGA_DUMP_LCP")[0] == '1');
    dump_records_as_tuples(out, &recv, with_lcp);
    fclose(out);

    record_buffer_free(&recv);
    Close_GDB(gdb);

    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == ROOT)
        fprintf(stderr, "[stage2] all ranks finished, outputs at %s.<rank>.tuples\n", argv[2]);

    MPI_Finalize();
    return 0;
}
