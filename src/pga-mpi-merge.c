// pga-mpi-merge: Phase 4b — distributed merge + alignment.
//
// Each rank runs:
//   1. Read_GDB, extract its share of records.
//   2. run_stage2 (Alltoallv redistribute by 10-bit prefix bucket + sort).
//   3. Build KmerStreamAdapter on local sorted slice.
//   4. Compute deterministic contig-pair → rank assignment (redundantly).
//   5. self_adaptamer_merge — emits seed pairs into per-destination-rank
//      files (N_Units[r], C_Units[r] for r in [0, world_size)).
//   6. MPI Alltoallv exchange of seed bytes by destination (N + C separately).
//   7. pair_sort_search on this rank's received seeds.
//   8. la_merge → ALNtoPAF.
//
// With world_size=1 the routing collapses to "everyone routes to rank 0,"
// which IS rank 0, and the output must match FastGA's gold output.
//
// Usage:  mpirun -n N pga-mpi-merge <gdb-stem>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mpi.h>

#include "extract.h"
#include "record.h"
#include "sort.h"
#include "dsort.h"
#include "run_cmd.h"
#include "constants.h"
#include "GDB.h"
#include "kmer_adapter.h"
#include "libfastk.h"
#include "align.h"
#include "alncode.h"
#include "gene_core.h"
#include "fastga_pipeline.h"
#include "contig_assignment.h"
#include "seed_exchange.h"
#include "chunked_exchange.h"

// constants.h's #define KMER 40 collides with the runtime int variable.
#undef KMER

#define ROOT 0

// Externs consumed by fastga_pipeline.c.
int    NTHREADS   = 1;
int    NPARTS     = 1;
int    SELF       = 1;
int    SOFT_MASK  = 0;
int    VERBOSE    = 1;
int    KMER       = 40;

int    FREQ        = 10;
int    CHAIN_BREAK = 2000;
int    CHAIN_MIN   = 170;
int    ALIGN_MIN   = 100;
double ALIGN_RATE  = 0.3;

int    IBYTE, JBYTE;
int    ICONT, JCONT;
int    IPOST, JPOST;
int    ISIGN, JSIGN;
int    KBYTE, CBYTE, LBYTE, PAYOFF;
int    ESHIFT;

int    NCONTS;
int64  AMXPOS, BMXPOS, MAXDAG;
int    DBYTE;

int    *Select   = NULL;
int    *IDBsplit = NULL;
int    *Perm1    = NULL;
int    *Perm2    = NULL;
char   *SORT_PATH = NULL;
char   *PAIR_NAME = NULL;
char   *ALGN_UNIQ = NULL;
char   *ALGN_PAIR = NULL;

IOBuffer *N_Units = NULL;
IOBuffer *C_Units = NULL;

SeedUnit *N_SeedUnits = NULL;
SeedUnit *C_SeedUnits = NULL;
uint8    *Recv_N_Data = NULL;
uint8    *Recv_C_Data = NULL;
int64     Recv_N_Len  = 0;
int64     Recv_C_Len  = 0;

// Phase 4b routing globals (consumed by fastga_pipeline.c).
int8_t *RankForPair = NULL;
int     NUNITS      = 1;

// --- Instrumentation (PR 1) ----------------------------------------------
//
// Memory: each stage's peak RSS comes from VmHWM in /proc/self/status. The
// high-water mark is reset at every stage boundary by writing "5" to
// /proc/self/clear_refs (Linux >= 4.0); if that fails the peaks are
// cumulative and reset_ok is 0 in the output.
//
// The reset also resets what getrusage / wait4 report as ru_maxrss, so
// /usr/bin/time -v and sacct MaxRSS then show only the last stage. Use
// overall_peak_kb from the [pr1-mem] line instead, or set
// PGA_MEM_STAGE_RESET=0 to keep the external numbers valid (per-stage peaks
// are then cumulative).

enum { MEM_EXTRACT, MEM_STAGE2, MEM_ADAPTER, MEM_MERGE, MEM_REDIST,
       MEM_ALIGN, MEM_CONVERT, MEM_NSTAGES };
static const char *Mem_Stage_Name[MEM_NSTAGES] =
  { "extract", "stage2", "adapter", "merge", "redist", "align", "convert" };
static long Mem_Peak_KB[MEM_NSTAGES];
static long Mem_End_KB[MEM_NSTAGES];
static int  Mem_Reset_OK = 1;

static void read_rss_kb(long *rss_kb, long *hwm_kb)
{
    char line[256];
    FILE *f = fopen("/proc/self/status", "r");
    *rss_kb = *hwm_kb = -1;
    if (f == NULL) return;
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "VmRSS:", 6) == 0) *rss_kb = atol(line + 6);
        else if (strncmp(line, "VmHWM:", 6) == 0) *hwm_kb = atol(line + 6);
    }
    fclose(f);
}

static void reset_hwm(void)
{
    const char *e = getenv("PGA_MEM_STAGE_RESET");
    if (e != NULL && strcmp(e, "0") == 0) { Mem_Reset_OK = 0; return; }
    int fd = open("/proc/self/clear_refs", O_WRONLY);
    if (fd < 0 || write(fd, "5", 1) != 1) Mem_Reset_OK = 0;
    if (fd >= 0) close(fd);
}

// Record the stage that just ended, then start the next one.
static void mem_stage_end(int stage)
{
    long rss, hwm;
    read_rss_kb(&rss, &hwm);
    Mem_Peak_KB[stage] = hwm;
    Mem_End_KB[stage]  = rss;
    reset_hwm();
}

// Per-task alignment records: one line per align_contigs call.
static FILE *Pair_Task_File = NULL;

static void write_pair_task(const PairTaskRecord *r)
{
    fprintf(Pair_Task_File,
            "%d\t%c\t%d\t%d\t%lld\t%.6f\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n",
            r->tid, r->comp ? 'C' : 'N', r->anchor, r->partner,
            (long long) r->nseeds, r->seconds, (long long) r->nhits,
            (long long) r->nlass, (long long) r->nlive, (long long) r->nlcov,
            (long long) r->ndiag, (long long) r->maxdiag,
            (long long) r->alen, (long long) r->blen);
}

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// Per-rank scratch directory for alignment files (.las, .1aln). Defaults to
// /tmp, which is memory-backed (tmpfs) on Perlmutter compute nodes; set
// PGA_SCRATCH_DIR (e.g. a $SCRATCH subdirectory) when those files would not
// fit in node memory alongside the process.
static char *make_temp_dir_for_rank(int rank)
{
    char tmpl[1024];
    const char *base = getenv("PGA_SCRATCH_DIR");
    if (base == NULL || base[0] == '\0') base = "/tmp";
    // Scratch file paths are built into 1024-byte buffers further on; refuse a
    // base long enough to be truncated anywhere rather than truncate silently.
    if (strlen(base) > 512) {
        fprintf(stderr, "rank %d: PGA_SCRATCH_DIR is longer than 512 characters\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    snprintf(tmpl, sizeof(tmpl), "%s/pga-mpi-merge-r%d-XXXXXX", base, rank);
    char *dir = mkdtemp(tmpl);
    if (dir == NULL) { perror("mkdtemp"); MPI_Abort(MPI_COMM_WORLD, 1); }
    return strdup(dir);
}

// Fill callback for chunked_alltoallv: copy from per-dest SeedUnit buffers.
typedef struct { SeedUnit *units; } SeedMemSource;

static void seed_mem_fill(void *ctx, int peer, int64_t off, int64_t len, uint8_t *dst)
{
    SeedMemSource *s = (SeedMemSource *) ctx;
    memcpy(dst, s->units[peer].data + off, len);
}

// In-memory seed exchange for one type (N or C).
// Returns a freshly malloc'd buffer with all received seeds in source-rank
// order (byte-identical to the old file-based exchange). Frees the send
// units before returning (memory rule 2).
static uint8 *redistribute_seeds_mem(SeedUnit *units, int world_size,
                                     int seed_bytes,
                                     int64 *recv_total_out,
                                     int64 *sent_out, XchgStats *xst_out)
{
    // 1-rank shortcut: no exchange needed, take ownership of the buffer.
    if (world_size == 1) {
        memset(xst_out, 0, sizeof(*xst_out));
        *sent_out = units[0].len;
        *recv_total_out = units[0].len;
        uint8 *buf = units[0].data;
        units[0].data = NULL;
        units[0].len = units[0].cap = 0;
        return buf;
    }

    int64_t *send_sizes  = calloc(world_size, sizeof(int64_t));
    int64_t *recv_sizes  = calloc(world_size, sizeof(int64_t));
    if (!send_sizes || !recv_sizes) {
        fprintf(stderr, "redistribute_seeds_mem: out of memory\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int64 send_total = 0;
    for (int r = 0; r < world_size; r++) {
        send_sizes[r] = units[r].len;
        send_total += send_sizes[r];
    }
    *sent_out = send_total;

    xchg_exchange_sizes(send_sizes, recv_sizes, MPI_COMM_WORLD);
    int64 recv_total = 0;
    for (int r = 0; r < world_size; r++)
        recv_total += recv_sizes[r];
    *recv_total_out = recv_total;

    uint8 *recv_buf = malloc(recv_total > 0 ? (size_t) recv_total : 1);
    if (recv_buf == NULL) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        fprintf(stderr, "[rank %d] redistribute_seeds_mem: malloc(%lld) failed\n",
                rank, (long long) recv_total);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    SeedMemSource sctx = { units };
    XchgSource xsrc = { NULL, seed_mem_fill, &sctx };
    XchgSink   xdst = { recv_buf, NULL, NULL };
    chunked_alltoallv(send_sizes, recv_sizes, seed_bytes,
                      xchg_round_budget_from_env(),
                      &xsrc, &xdst, MPI_COMM_WORLD, xst_out);

    for (int r = 0; r < world_size; r++)
        seed_unit_free(&units[r]);

    free(send_sizes);
    free(recv_sizes);
    return recv_buf;
}

// Count per-A-contig seeds in a received buffer (replaces recompute_buck).
static void count_buck(int64 *buck, const uint8 *buf, int64 total, int seed_bytes)
{
    int   icont_off = 1 + IPOST;
    int64 mask      = ((int64) 1 << (8 * ICONT - 1)) - 1;
    bzero(buck, sizeof(int64) * NCONTS);
    for (int64 i = 0; i < total; i += seed_bytes) {
        int64 icont = 0;
        memcpy(&icont, buf + i + icont_off, ICONT);
        icont &= mask;
        buck[icont]++;
    }
}

int main(int argc, char *argv[])
{
    Prog_Name = "pga-mpi-merge";

    MPI_Init(&argc, &argv);
    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc != 2) {
        if (my_rank == ROOT)
            fprintf(stderr, "Usage: %s <gdb-stem>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }
    char *gdb_stem = argv[1];

    GDB _gdb, *gdb = &_gdb;
    if (Read_GDB(gdb, gdb_stem) < 0) {
        fprintf(stderr, "rank %d: Read_GDB failed for '%s'\n", my_rank, gdb_stem);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (my_rank == ROOT) {
        fprintf(stderr, "[merge] world_size=%d  GDB: %d contigs, %lld total bases\n",
                world_size, gdb->ncontig, (long long) gdb->seqtot);
    }

    RecordSizing sizing = compute_record_sizing(gdb);

    // Stage timers (max reduced and printed at the end).
    double t_extract = 0, t_sort = 0, t_merge = 0, t_alignments = 0;
    double t_stage_start;

    // ---- STAGE: extract -------------------------------------------------
    reset_hwm();
    t_stage_start = MPI_Wtime();
    RecordBuffer local;
    record_buffer_init(&local, sizing);
    extract_stats_reset();
    extract_to_records(gdb, world_size, my_rank, &local);
    t_extract = MPI_Wtime() - t_stage_start;
    int64_t extracted_records = local.count;

    // ---- STAGE: sort ----------------------------------------------------
    t_stage_start = MPI_Wtime();
    RecordBuffer recv;
    memset(&recv, 0, sizeof(recv));
    mem_stage_end(MEM_EXTRACT);
    run_stage2(&local, &recv, MPI_COMM_WORLD);
    record_buffer_free(&local);
    t_sort = MPI_Wtime() - t_stage_start;
    mem_stage_end(MEM_STAGE2);

    fprintf(stderr, "[merge] rank %d: own slice has %lld records (sorted)\n",
            my_rank, (long long) recv.count);

    // Per-rank extraction/stage-2 breakdown (O1-0). One line per rank,
    // key=value, so scripts can grep '^\[o1-0\]'.
    {
        const ExtractStats *es = &g_extract_stats;
        const Stage2Timing *st = &g_stage2_timing;
        fprintf(stderr,
                "[o1-0] rank %d: extract_s=%.6f contigs=%lld bases=%lld "
                "records=%lld window_reads=%lld window_bases=%lld "
                "pack_reads=%lld pack_bases=%lld window_packs=%lld "
                "hist_s=%.6f pack_s=%.6f exchange_s=%.6f shuffle_s=%.6f "
                "local_sort_s=%.6f verify_s=%.6f sort_stage_s=%.6f "
                "slice_records=%lld\n",
                my_rank, t_extract, (long long) es->contigs,
                (long long) es->bases, (long long) extracted_records,
                (long long) es->window_reads, (long long) es->window_bases,
                (long long) es->pack_reads, (long long) es->pack_bases,
                (long long) es->window_packs,
                st->hist_s, st->pack_s, st->exchange_s, st->shuffle_s,
                st->local_sort_s, st->verify_s, t_sort,
                (long long) recv.count);
    }

    // ---- STAGE: merge ---------------------------------------------------
    // Covers: build_kmer_adapter + extern setup + IOBuffer/SeedUnit alloc +
    // self_adaptamer_merge + redistribute_seeds_mem(N+C) + count_buck.
    t_stage_start = MPI_Wtime();

    KmerStreamAdapter T1;
    PostListAdapter   P1;
    memset(&T1, 0, sizeof(T1));
    memset(&P1, 0, sizeof(P1));
    if (build_kmer_adapter(&recv, &T1, &P1) < 0) {
        fprintf(stderr, "rank %d: build_kmer_adapter failed\n", my_rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    fprintf(stderr, "[merge] rank %d: adapter %lld entries, pbyte=%d\n",
            my_rank, (long long) T1.nels, T1.pbyte);
    double t_adapter = MPI_Wtime() - t_stage_start;
    int64 slice_records = recv.count;
    int64 record_bytes  = recv.count * (int64) recv.sizing.record_size;
    int64 adapter_bytes = T1.nels * (int64) T1.pbyte + (int64) T1.ixlen * (int64) sizeof(int64);

    // The adapter holds its own copy of every record (build_kmer_adapter
    // copies, it does not alias), and nothing after this point reads recv.
    record_buffer_free(&recv);
    mem_stage_end(MEM_ADAPTER);

    // Set FastGA externs.
    KMER   = T1.kmer;
    IBYTE  = P1.pbyte;
    ICONT  = P1.cbyte;
    IPOST  = IBYTE - ICONT;
    ISIGN  = IBYTE - 1;
    JBYTE  = IBYTE;
    JCONT  = ICONT;
    JPOST  = IPOST;
    JSIGN  = ISIGN;
    KBYTE  = T1.pbyte;
    CBYTE  = T1.hbyte;
    LBYTE  = CBYTE + 1;
    PAYOFF = LBYTE + 1;
    ESHIFT = 8 * IPOST;
    NCONTS = gdb->ncontig;

    AMXPOS = 0;
    for (int r = 0; r < gdb->ncontig; r++) {
        int64 len = gdb->contigs[r].clen;
        if (len > AMXPOS) AMXPOS = len;
    }
    BMXPOS = AMXPOS;
    MAXDAG = AMXPOS + BMXPOS;

    DBYTE = 0;
    int64 cum = 1;
    while (cum < MAXDAG) { cum *= 256; DBYTE++; }

    Perm1 = malloc(gdb->ncontig * sizeof(int));
    Perm2 = Perm1;
    for (int i = 0; i < gdb->ncontig; i++) Perm1[i] = i;

    NPARTS = 1;
    IDBsplit = malloc((NPARTS + 1) * sizeof(int));
    IDBsplit[0] = 0;
    IDBsplit[NPARTS] = gdb->ncontig;

    Select = malloc(gdb->ncontig * sizeof(int));
    for (int i = 0; i < gdb->ncontig; i++) Select[i] = 0;

    ContigAssignment assign;
    if (contig_assignment_build(gdb, world_size, &assign) < 0) {
        fprintf(stderr, "rank %d: contig_assignment_build failed\n", my_rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (my_rank == ROOT) contig_assignment_print(&assign, gdb);

    // Make the assignment visible to fastga_pipeline.c routing.
    // We allocate a flat NCONTS*NCONTS map of int8_t, with both (i,j) and
    // (j,i) populated so the merge code can read RankForPair[icont*NCONTS+jcont]
    // without a min/max swap on every seed.
    RankForPair = malloc((size_t) NCONTS * NCONTS);
    if (RankForPair == NULL) {
        fprintf(stderr, "rank %d: out of memory for RankForPair\n", my_rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (int i = 0; i < NCONTS; i++) {
        for (int j = 0; j < NCONTS; j++) {
            int lo = i < j ? i : j;
            int hi = i < j ? j : i;
            RankForPair[(int64) i * NCONTS + j] =
                assign.rank_for_pair[(int64) lo * NCONTS + hi];
        }
    }
    NUNITS = world_size;

    SORT_PATH = make_temp_dir_for_rank(my_rank);
    fprintf(stderr, "[merge] rank %d: scratch %s\n", my_rank, SORT_PATH);

    PAIR_NAME = Malloc(256, "pair_name");
    ALGN_UNIQ = Malloc(256, "algn_uniq");
    ALGN_PAIR = Malloc(256, "algn_pair");
    sprintf(PAIR_NAME, "_pair.%d", getpid());
    sprintf(ALGN_UNIQ, "_uniq.%d", getpid());
    sprintf(ALGN_PAIR, "_algn.%d", getpid());

    // NUNITS = world_size slots per local thread, matching the per-destination
    // routing in new_self_merge_thread.
    int num_units = NTHREADS * NUNITS;
    N_Units = calloc(num_units, sizeof(IOBuffer));
    C_Units = calloc(num_units, sizeof(IOBuffer));
    int64 *unit_bucks  = calloc(2ll * NTHREADS * NCONTS, sizeof(int64));
    if (N_Units == NULL || C_Units == NULL || unit_bucks == NULL) {
        fprintf(stderr, "rank %d: out of memory for unit buffers\n", my_rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    N_SeedUnits = calloc(num_units, sizeof(SeedUnit));
    C_SeedUnits = calloc(num_units, sizeof(SeedUnit));
    if (N_SeedUnits == NULL || C_SeedUnits == NULL) {
        fprintf(stderr, "rank %d: out of memory for SeedUnits\n", my_rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int k = 0;
    for (int i = 0; i < NTHREADS; i++) {
        for (int j = 0; j < NUNITS; j++) {
            N_Units[k].buck = unit_bucks + (2ll * i) * NCONTS;
            C_Units[k].buck = unit_bucks + (2ll * i + 1) * NCONTS;
            N_Units[k].inum = k;
            C_Units[k].inum = k;
            seed_unit_init(&N_SeedUnits[k], k);
            seed_unit_init(&C_SeedUnits[k], k);
            k++;
        }
    }

    if (my_rank == ROOT)
        fprintf(stderr, "[merge] starting self-adaptamer merge on each rank...\n");

    // self_adaptamer_merge expects gdb->seqtot but we want the per-rank
    // record count proxy here for the verbose output ratio. Pass seqtot
    // since it's only used for a logging division.
    double t0 = now_s();
    self_adaptamer_merge((Kmer_Stream *)&T1, (Post_List *)&P1, gdb->seqtot);
    double t_merge_core = now_s() - t0;

    // Seeds are now in SeedUnits; pair_sort_search and everything after it
    // read only those buffers and the GDB, never the k-mer adapter.
    free_kmer_adapter(&T1, &P1);
    mem_stage_end(MEM_MERGE);

    // Trim slack from doubling before exchange.
    for (int r = 0; r < num_units; r++) {
        seed_unit_trim(&N_SeedUnits[r]);
        seed_unit_trim(&C_SeedUnits[r]);
    }

    if (my_rank == ROOT)
        fprintf(stderr, "[merge] redistributing seeds via MPI_Alltoallv...\n");

    int seed_bytes = 1 + IBYTE + JBYTE;
    int64 n_sent, c_sent, n_recv, c_recv;
    XchgStats n_xst, c_xst;
    t0 = now_s();
    Recv_N_Data = redistribute_seeds_mem(N_SeedUnits, world_size, seed_bytes,
                                         &n_recv, &n_sent, &n_xst);
    Recv_N_Len = n_recv;
    double t_redist_n = now_s() - t0;
    t0 = now_s();
    Recv_C_Data = redistribute_seeds_mem(C_SeedUnits, world_size, seed_bytes,
                                         &c_recv, &c_sent, &c_xst);
    Recv_C_Len = c_recv;
    double t_redist_c = now_s() - t0;
    fprintf(stderr, "[merge] rank %d: received N=%lld C=%lld bytes\n",
            my_rank, (long long) n_recv, (long long) c_recv);

    // After exchange we have ONE slot per type. NUNITS for downstream code
    // is now 1 (matching pair_sort_search's NPARTS=1 expectation).
    NUNITS = 1;

    // Recount per-A-contig seeds from received buffers; on 1 rank the merge
    // counts are already correct (all seeds stayed local).
    t0 = now_s();
    if (world_size > 1) {
        count_buck(N_Units[0].buck, Recv_N_Data, n_recv, seed_bytes);
        count_buck(C_Units[0].buck, Recv_C_Data, c_recv, seed_bytes);
    }
    double t_recompute = now_s() - t0;
    t_merge = MPI_Wtime() - t_stage_start;
    mem_stage_end(MEM_REDIST);

    // ---- STAGE: alignments ----------------------------------------------
    // Covers pair_sort_search + ALNtoPAF.
    t_stage_start = MPI_Wtime();

    // pair_sort_search uses NPARTS=1, so it only looks at N_Units[0] and
    // C_Units[0] — exactly the slots we filled with received bytes.
    fprintf(stderr, "[merge] rank %d: starting pair_sort_search\n", my_rank);
    // Per-task records (pairs.rank<r>.tsv in the working directory) only when
    // PGA_PAIR_TASKS is set to something other than "0"; off by default.
    const char *pt_env = getenv("PGA_PAIR_TASKS");
    if (pt_env != NULL && pt_env[0] != '\0' && strcmp(pt_env, "0") != 0) {
        char task_path[64];
        snprintf(task_path, sizeof(task_path), "pairs.rank%d.tsv", my_rank);
        Pair_Task_File = fopen(task_path, "w");
        if (Pair_Task_File != NULL) {
            fprintf(Pair_Task_File, "tid\ttype\tanchor\tpartner\tseeds\tseconds\tnhits\tnlass"
                                    "\tnlive\tnlcov\tndiag\tmaxdiag\talen\tblen\n");
            Pair_Task_Hook = write_pair_task;
        } else
            fprintf(stderr, "rank %d: cannot open %s; per-task records disabled\n",
                    my_rank, task_path);
    }
    t0 = now_s();
    pair_sort_search(gdb, gdb);
    double t_pair_search = now_s() - t0;
    Pair_Task_Hook = NULL;
    if (Pair_Task_File != NULL) fclose(Pair_Task_File);
    mem_stage_end(MEM_ALIGN);

    // Convert with a time limit; keep ALNtoPAF's stderr. A failure on any
    // rank makes the run exit non-zero (see the end of main).
    char aln_path[1024], paf_path[64], err_path[64];
    snprintf(aln_path, sizeof(aln_path), "%s/%s.1aln", SORT_PATH, ALGN_UNIQ);
    snprintf(paf_path, sizeof(paf_path), "out.rank%d.paf", my_rank);
    snprintf(err_path, sizeof(err_path), "out.rank%d.alntopaf.err", my_rank);
    double conv_limit = alntopaf_timeout_from_env(600.0);
    char *conv_argv[] = { "ALNtoPAF", aln_path, NULL };
    fprintf(stderr, "[merge] rank %d: ALNtoPAF %s > %s 2> %s (limit %.0f s)\n",
            my_rank, aln_path, paf_path, err_path, conv_limit);
    RunResult conv = run_with_timeout(conv_argv, paf_path, err_path, conv_limit, 10.0);
    int conv_failed = (conv.status != RUN_OK);
    if (conv_failed)
    {
        char why[256], partial[80];
        run_result_describe(&conv, conv_limit, why, sizeof(why));
        fprintf(stderr, "rank %d: ALNtoPAF failed: %s; stderr in %s\n",
                my_rank, why, err_path);
        snprintf(partial, sizeof(partial), "%s.partial", paf_path);
        rename(paf_path, partial);
    }
    else
    {
        struct stat st;
        if (stat(err_path, &st) == 0 && st.st_size == 0)
            unlink(err_path);
    }
    t_alignments = MPI_Wtime() - t_stage_start;
    mem_stage_end(MEM_CONVERT);

    // Per-rank merge/alignment breakdown (PR 1). Seeds are 1+IBYTE+JBYTE bytes.
    {
        const PairSearchStats *ps = &Pair_Search_Stats;
        int seed_bytes = 1 + IBYTE + JBYTE;
        fprintf(stderr,
                "[pr1] rank %d: slice_records=%lld adapter_s=%.6f merge_core_s=%.6f "
                "redist_n_s=%.6f redist_c_s=%.6f recompute_s=%.6f "
                "pair_search_s=%.6f reimport_n_s=%.6f reimport_c_s=%.6f "
                "sort_n_s=%.6f sort_c_s=%.6f search_n_s=%.6f search_c_s=%.6f "
                "la_sort_s=%.6f la_merge_s=%.6f convert_s=%.6f "
                "seed_bytes=%d sent_n_bytes=%lld sent_c_bytes=%lld "
                "recv_n_bytes=%lld recv_c_bytes=%lld seeds_n=%lld seeds_c=%lld "
                "record_bytes=%lld adapter_bytes=%lld sarray_bytes=%lld "
                "stage2_rounds=%lld stage2_staging=%lld redist_n_rounds=%lld "
                "redist_c_rounds=%lld redist_staging=%lld\n",
                my_rank, (long long) slice_records, t_adapter, t_merge_core,
                t_redist_n, t_redist_c, t_recompute,
                t_pair_search, ps->reimport_s[0], ps->reimport_s[1],
                ps->sort_s[0], ps->sort_s[1], ps->search_s[0], ps->search_s[1],
                ps->la_sort_s, ps->la_merge_s, conv.elapsed_s,
                seed_bytes, (long long) n_sent, (long long) c_sent,
                (long long) n_recv, (long long) c_recv,
                (long long) ps->nels[0], (long long) ps->nels[1],
                (long long) record_bytes, (long long) adapter_bytes,
                (long long) ps->sarray_bytes,
                (long long) g_stage2_timing.exchange_rounds,
                (long long) g_stage2_timing.exchange_staging,
                (long long) n_xst.rounds, (long long) c_xst.rounds,
                (long long) (n_xst.staging_bytes > c_xst.staging_bytes
                             ? n_xst.staging_bytes : c_xst.staging_bytes));
        char memline[512];
        long overall = 0;
        for (int s = 0; s < MEM_NSTAGES; s++)
            if (Mem_Peak_KB[s] > overall) overall = Mem_Peak_KB[s];
        int off = snprintf(memline, sizeof(memline),
                           "[pr1-mem] rank %d: reset_ok=%d overall_peak_kb=%ld",
                           my_rank, Mem_Reset_OK, overall);
        for (int s = 0; s < MEM_NSTAGES && off < (int) sizeof(memline); s++)
            off += snprintf(memline + off, sizeof(memline) - off,
                            " %s_peak_kb=%ld %s_end_kb=%ld",
                            Mem_Stage_Name[s], Mem_Peak_KB[s],
                            Mem_Stage_Name[s], Mem_End_KB[s]);
        fprintf(stderr, "%s\n", memline);
    }

    // Reduce stage timings to max-across-ranks; print + CSV on rank 0.
    {
        double local_t[4] = { t_extract, t_sort, t_merge, t_alignments };
        double max_t[4]   = { 0, 0, 0, 0 };
        MPI_Reduce(local_t, max_t, 4, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
        if (my_rank == ROOT) {
            double total = max_t[0] + max_t[1] + max_t[2] + max_t[3];
            fprintf(stderr,
                    "\n=== Stage timing (max across %d ranks) ===\n"
                    "  extract     %8.3fs\n"
                    "  sort        %8.3fs\n"
                    "  merge       %8.3fs\n"
                    "  alignments  %8.3fs\n"
                    "  TOTAL       %8.3fs\n",
                    world_size,
                    max_t[0], max_t[1], max_t[2], max_t[3], total);

            // Append CSV row for the test driver to consume.
            const char *csv_path = "pga-mpi-merge.timing.csv";
            FILE *csv_f = fopen(csv_path, "a");
            if (csv_f != NULL) {
                if (ftell(csv_f) == 0) {
                    fprintf(csv_f, "world_size,extract_s,sort_s,merge_s,"
                                   "alignments_s,total_s\n");
                }
                fprintf(csv_f, "%d,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                        world_size, max_t[0], max_t[1], max_t[2],
                        max_t[3], total);
                fclose(csv_f);
            }
        }
    }

cleanup:
    free_kmer_adapter(&T1, &P1);
    record_buffer_free(&recv);
    contig_assignment_free(&assign);
    free(RankForPair);
    Close_GDB(gdb);

    // Remove the scratch directory without a shell: the path may contain
    // spaces or shell metacharacters now that PGA_SCRATCH_DIR sets its base.
    char *rm_argv[] = { "rm", "-rf", "--", SORT_PATH, NULL };
    RunResult rm_res = run_with_timeout(rm_argv, NULL, NULL, 300.0, 10.0);
    if (rm_res.status != RUN_OK)
        fprintf(stderr, "rank %d: could not remove scratch directory %s\n",
                my_rank, SORT_PATH);

    free(SORT_PATH);
    free(PAIR_NAME);
    free(ALGN_UNIQ);
    free(ALGN_PAIR);
    free(Perm1);
    free(IDBsplit);
    free(Select);
    free(unit_bucks);
    free(N_Units);
    free(C_Units);
    free(N_SeedUnits);
    free(C_SeedUnits);

    int any_conv_failed = 0;
    MPI_Allreduce(&conv_failed, &any_conv_failed, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == ROOT)
    {
        if (any_conv_failed)
            fprintf(stderr, "[merge] FAILED: ALNtoPAF failed on at least one "
                            "rank; PAF output is incomplete\n");
        else
            fprintf(stderr, "[merge] done\n");
    }
    MPI_Finalize();
    return any_conv_failed ? 4 : 0;
}
