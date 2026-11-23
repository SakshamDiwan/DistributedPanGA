// pga-dsort-bench: extraction + distributed sort only (R3 experiments).
//
// Runs exactly the production path of pga-mpi-merge up to the end of the sort
// stage (extract_to_records, barrier, run_stage2 / run_stage2_take), with a
// barrier on both sides of the sort so the stage's wall time is clean, then
// prints per-rank lines and three global digests and exits. No alignment.
//
//   [o1-0] / [r3-s2] rank lines   as in pga-mpi-merge
//   [bench] rep=k ...             rank 0: wall times, totals, digests
//
// Digests (64-bit sums over all records, independent of rank count):
//   kseq   sum_i H(i, k-mer bytes, LCP byte), i = global sorted index. Must be
//          identical for every rank count and every sort variant: it fixes the
//          sorted k-mer sequence and the semantic LCP of every record.
//   multi  sum H(record without LCP byte): the record multiset (nothing lost
//          or duplicated).
//   order  sum_i H(i, record without LCP byte): also fixes the order of records
//          with equal k-mers. Informational: that order is not specified by
//          msd_sort and already differs between rank counts.
//
// Usage: srun ... pga-dsort-bench <gdb-stem> [reps]
// Env: PGA_STAGE2_KEEP_LOCAL=1 (old buffer lifetime), PGA_SORT_THREADS,
//      PGA_SORT_PARTS, PGA_EXCHANGE_ROUND_BYTES, PGA_BUCKET_HIST.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>
#include <mpi.h>

#include "extract.h"
#include "record.h"
#include "sort.h"
#include "dsort.h"
#include "constants.h"
#include "GDB.h"

static uint64_t mix64(uint64_t z)
{
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static uint64_t hash_bytes(const uint8_t *p, int n, uint64_t seed)
{
    uint64_t h = 0xcbf29ce484222325ULL ^ seed;
    for (int i = 0; i < n; i++) { h ^= p[i]; h *= 0x100000001b3ULL; }
    return mix64(h);
}

static void read_mem(long *rss, long *hwm)
{
    char line[256];
    FILE *f = fopen("/proc/self/status", "r");
    *rss = *hwm = -1;
    if (f == NULL) return;
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "VmRSS:", 6) == 0) *rss = atol(line + 6);
        else if (strncmp(line, "VmHWM:", 6) == 0) *hwm = atol(line + 6);
    }
    fclose(f);
}

static void reset_hwm(void)
{
    int fd = open("/proc/self/clear_refs", O_WRONLY);
    if (fd >= 0) { if (write(fd, "5", 1) != 1) { /* peaks stay cumulative */ } close(fd); }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc != 2 && argc != 3) {
        if (my_rank == 0) fprintf(stderr, "Usage: %s <gdb-stem> [reps]\n", argv[0]);
        MPI_Finalize();
        return 1;
    }
    int reps = (argc == 3) ? atoi(argv[2]) : 1;
    if (reps < 1) reps = 1;

    GDB _gdb, *gdb = &_gdb;
    if (Read_GDB(gdb, argv[1]) != 0) {
        fprintf(stderr, "rank %d: Read_GDB failed for '%s'\n", my_rank, argv[1]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    RecordSizing sizing = compute_record_sizing(gdb);
    const char *keep_env = getenv("PGA_STAGE2_KEEP_LOCAL");
    int keep = (keep_env != NULL && strcmp(keep_env, "1") == 0);
    if (my_rank == 0)
        fprintf(stderr, "[bench] world=%d contigs=%d bases=%lld rsize=%d keep_local=%d "
                        "sort_threads=%d sort_parts=%d\n",
                world_size, gdb->ncontig, (long long) gdb->seqtot, sizing.record_size,
                keep, sort_threads_from_env(), sort_parts_from_env());

    for (int rep = 1; rep <= reps; rep++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        reset_hwm();
        double t0 = MPI_Wtime();
        RecordBuffer local;
        record_buffer_init(&local, sizing);
        extract_stats_reset();
        extract_to_records(gdb, world_size, my_rank, &local);
        double t_ext_own = MPI_Wtime() - t0;
        int64_t extracted = local.count;
        int64_t cap_bytes = local.capacity * (int64_t) sizing.record_size;
        long ext_rss, ext_hwm;
        read_mem(&ext_rss, &ext_hwm);
        MPI_Barrier(MPI_COMM_WORLD);
        double t1 = MPI_Wtime();

        reset_hwm();
        RecordBuffer recv;
        memset(&recv, 0, sizeof(recv));
        if (keep)
            run_stage2(&local, &recv, MPI_COMM_WORLD);
        else
            run_stage2_take(&local, &recv, MPI_COMM_WORLD);
        double t_sort_own = MPI_Wtime() - t1;
        record_buffer_free(&local);
        long s2_rss, s2_hwm;
        read_mem(&s2_rss, &s2_hwm);
        MPI_Barrier(MPI_COMM_WORLD);
        double t2 = MPI_Wtime();

        const ExtractStats *es = &g_extract_stats;
        const Stage2Timing *st = &g_stage2_timing;
        fprintf(stderr,
                "[o1-0] rank %d: extract_s=%.6f contigs=%lld bases=%lld records=%lld "
                "window_reads=%lld window_bases=%lld pack_reads=%lld pack_bases=%lld "
                "window_packs=%lld hist_s=%.6f pack_s=%.6f exchange_s=%.6f shuffle_s=%.6f "
                "local_sort_s=%.6f verify_s=%.6f sort_stage_s=%.6f slice_records=%lld "
                "rep=%d ext_cap_bytes=%lld ext_rss_kb=%ld ext_hwm_kb=%ld s2_hwm_kb=%ld\n",
                my_rank, t_ext_own, (long long) es->contigs, (long long) es->bases,
                (long long) extracted, (long long) es->window_reads,
                (long long) es->window_bases, (long long) es->pack_reads,
                (long long) es->pack_bases, (long long) es->window_packs,
                st->hist_s, st->pack_s, st->exchange_s, st->shuffle_s, st->local_sort_s,
                st->verify_s, t_sort_own, (long long) recv.count, rep,
                (long long) cap_bytes, ext_rss, ext_hwm, s2_hwm);
        stage2_report(stderr, my_rank);

        // Digests.
        int64_t first = 0;
        MPI_Exscan(&recv.count, &first, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
        if (my_rank == 0) first = 0;
        int rs = sizing.record_size;
        uint64_t d[3] = { 0, 0, 0 };
        for (int64_t i = 0; i < recv.count; i++)
        {
            const uint8_t *r = recv.data + i * rs;
            uint64_t gi = (uint64_t) (first + i);
            d[0] += hash_bytes(r, 1 + KBYTES, mix64(gi + 1));   // LCP + k-mer
            d[1] += hash_bytes(r + 1, rs - 1, 0);                // record w/o LCP
            d[2] += hash_bytes(r + 1, rs - 1, mix64(gi + 1));
        }
        uint64_t g[3];
        MPI_Reduce(d, g, 3, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

        double loc[4] = { t_ext_own, t_sort_own, (double) s2_hwm, (double) ext_hwm }, mx[4];
        MPI_Reduce(loc, mx, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        double sum_hwm = 0, l_hwm = (double) (s2_hwm > ext_hwm ? s2_hwm : ext_hwm);
        MPI_Reduce(&l_hwm, &sum_hwm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        int64_t tot = 0;
        MPI_Reduce(&recv.count, &tot, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
        if (my_rank == 0)
            fprintf(stderr,
                    "[bench] rep=%d world=%d records=%lld extract_wall_s=%.6f sort_wall_s=%.6f "
                    "extract_max_own_s=%.6f sort_max_own_s=%.6f max_s2_hwm_kb=%.0f "
                    "max_ext_hwm_kb=%.0f sum_peak_kb=%.0f kseq=%016llx multi=%016llx order=%016llx\n",
                    rep, world_size, (long long) tot, t1 - t0, t2 - t1, mx[0], mx[1], mx[2], mx[3],
                    sum_hwm, (unsigned long long) g[0], (unsigned long long) g[1],
                    (unsigned long long) g[2]);
        record_buffer_free(&recv);
    }

    Close_GDB(gdb);
    MPI_Finalize();
    return 0;
}
