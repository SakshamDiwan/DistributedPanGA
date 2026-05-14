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
// which IS rank 0, and the output is byte-identical to single-rank pga-merge.
//
// Usage:  mpirun -n N pga-mpi-merge <gdb-stem>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>
#include <limits.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mpi.h>

#include "extract.h"
#include "record.h"
#include "sort.h"
#include "dsort.h"
#include "constants.h"
#include "GDB.h"
#include "kmer_adapter.h"
#include "libfastk.h"
#include "align.h"
#include "alncode.h"
#include "gene_core.h"
#include "fastga_pipeline.h"
#include "contig_assignment.h"

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

// Phase 4b routing globals (consumed by fastga_pipeline.c).
int8_t *RankForPair = NULL;
int     NUNITS      = 1;

static char *make_temp_dir_for_rank(int rank)
{
    char tmpl[64];
    snprintf(tmpl, sizeof(tmpl), "/tmp/pga-mpi-merge-r%d-XXXXXX", rank);
    char *dir = mkdtemp(tmpl);
    if (dir == NULL) { perror("mkdtemp"); MPI_Abort(MPI_COMM_WORLD, 1); }
    return strdup(dir);
}

// Read entire file into a freshly malloc'd buffer.
// Caller frees. Returns size in bytes; out_buf gets the buffer pointer.
static int64 read_file_all(int fd, uint8 **out_buf)
{
    off_t size = lseek(fd, 0, SEEK_END);
    if (size < 0) { perror("lseek"); MPI_Abort(MPI_COMM_WORLD, 1); }
    if (lseek(fd, 0, SEEK_SET) < 0) { perror("lseek"); MPI_Abort(MPI_COMM_WORLD, 1); }
    uint8 *buf = malloc(size > 0 ? (size_t) size : 1);
    if (buf == NULL) { perror("malloc"); MPI_Abort(MPI_COMM_WORLD, 1); }
    int64 total = 0;
    while (total < size) {
        ssize_t n = read(fd, buf + total, size - total);
        if (n <= 0) { perror("read"); MPI_Abort(MPI_COMM_WORLD, 1); }
        total += n;
    }
    *out_buf = buf;
    return size;
}

// Redistribute seeds for a single seed type (N or C) via MPI Alltoallv.
//
// Inputs:
//   units       — N_Units or C_Units, length world_size (one outgoing file per dest)
//   world_size  — # MPI ranks
//   tag         — 'N' or 'C', for diagnostics
//   sort_path   — scratch dir
//   pair_name   — scratch file basename
//   out_path    — caller-provided buffer; receives path of fresh local file
// Returns:
//   total bytes received (also = size of file at *out_path).
//
// Side effect: closes/replaces units[0].file with the combined incoming file.
static int64 redistribute_seeds(IOBuffer *units, int world_size, char tag,
                                const char *sort_path, const char *pair_name,
                                char *out_path, size_t out_path_len)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Single-rank shortcut: nothing to exchange. Just keep units[0].file
    // (which already holds all our seeds, since RankForPair routed every
    // pair to rank 0 = ourselves). Skip the Alltoallv to avoid the INT_MAX
    // chunking issue when total seed bytes > 2 GB.
    if (world_size == 1) {
        off_t size = lseek(units[0].file, 0, SEEK_END);
        if (lseek(units[0].file, 0, SEEK_SET) < 0) {
            perror("lseek"); MPI_Abort(MPI_COMM_WORLD, 1);
        }
        snprintf(out_path, out_path_len, "%s/%s.0.%c",
                 sort_path, pair_name, tag);
        return size;
    }

    int64 *send_sizes64 = calloc(world_size, sizeof(int64));
    int64  send_total64 = 0;
    for (int r = 0; r < world_size; r++) {
        send_sizes64[r] = lseek(units[r].file, 0, SEEK_END);
        if (send_sizes64[r] < 0) { perror("lseek"); MPI_Abort(MPI_COMM_WORLD, 1); }
        send_total64 += send_sizes64[r];
        if (send_sizes64[r] > (int64) INT_MAX) {
            fprintf(stderr, "rank %d: outgoing %c bytes to rank %d (%lld) > INT_MAX, "
                            "MPI_Alltoallv chunking required (TODO)\n",
                    my_rank, tag, r, (long long) send_sizes64[r]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    uint8 *send_buf = malloc(send_total64 > 0 ? (size_t) send_total64 : 1);
    if (send_buf == NULL) { perror("malloc send_buf"); MPI_Abort(MPI_COMM_WORLD, 1); }

    int *send_counts = calloc(world_size, sizeof(int));
    int *send_displs = calloc(world_size, sizeof(int));
    int64 cum = 0;
    for (int r = 0; r < world_size; r++) {
        send_displs[r] = (int) cum;
        send_counts[r] = (int) send_sizes64[r];
        if (send_sizes64[r] > 0) {
            if (lseek(units[r].file, 0, SEEK_SET) < 0) { perror("lseek"); MPI_Abort(MPI_COMM_WORLD, 1); }
            int64 got = 0;
            while (got < send_sizes64[r]) {
                ssize_t n = read(units[r].file, send_buf + cum + got, send_sizes64[r] - got);
                if (n <= 0) { perror("read"); MPI_Abort(MPI_COMM_WORLD, 1); }
                got += n;
            }
        }
        cum += send_sizes64[r];
        close(units[r].file);
        char path[1024];
        snprintf(path, sizeof(path), "%s/%s.%d.%c",
                 sort_path, pair_name, units[r].inum, tag);
        unlink(path);
    }

    int *recv_counts = calloc(world_size, sizeof(int));
    int *recv_displs = calloc(world_size, sizeof(int));
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);

    int64 recv_total64 = 0;
    for (int r = 0; r < world_size; r++) {
        recv_displs[r] = (int) recv_total64;
        recv_total64 += recv_counts[r];
        if (recv_total64 > (int64) INT_MAX) {
            fprintf(stderr, "rank %d: incoming %c total (%lld) > INT_MAX, "
                            "MPI_Alltoallv chunking required (TODO)\n",
                    my_rank, tag, (long long) recv_total64);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    uint8 *recv_buf = malloc(recv_total64 > 0 ? (size_t) recv_total64 : 1);
    if (recv_buf == NULL) { perror("malloc recv_buf"); MPI_Abort(MPI_COMM_WORLD, 1); }

    MPI_Alltoallv(send_buf, send_counts, send_displs, MPI_BYTE,
                  recv_buf, recv_counts, recv_displs, MPI_BYTE,
                  MPI_COMM_WORLD);

    free(send_buf);
    free(send_sizes64);
    free(send_counts);
    free(send_displs);
    free(recv_counts);
    free(recv_displs);

    // Write recv_buf to a fresh local file (will become units[0].file).
    snprintf(out_path, out_path_len, "%s/%s.combined.%c",
             sort_path, pair_name, tag);
    int fd = open(out_path, O_RDWR | O_CREAT | O_TRUNC, 0666);
    if (fd < 0) { perror("open combined"); MPI_Abort(MPI_COMM_WORLD, 1); }
    int64 written = 0;
    while (written < recv_total64) {
        ssize_t n = write(fd, recv_buf + written, recv_total64 - written);
        if (n <= 0) { perror("write combined"); MPI_Abort(MPI_COMM_WORLD, 1); }
        written += n;
    }
    free(recv_buf);

    units[0].file = fd;
    units[0].inum = 0;
    return recv_total64;
}

// Recompute units[0].buck[icont] = number of seeds (per A-contig) in the
// combined file at units[0].file. Needed because pair_sort_search consumes
// these counts and after redistribution they no longer match what merge
// originally wrote.
static void recompute_buck(IOBuffer *unit, int ibyte, int ipost, int icont_bytes)
{
    int fd = unit->file;
    off_t size = lseek(fd, 0, SEEK_END);
    if (lseek(fd, 0, SEEK_SET) < 0) { perror("lseek recompute"); MPI_Abort(MPI_COMM_WORLD, 1); }

    // Each seed: 1 byte LCP + IBYTE pay1 + JBYTE p
    // Within pay1: IPOST bytes pos + ICONT bytes contig.  The icont byte
    // is at offset 1 + IPOST inside the seed record.
    int seed_bytes = ibyte + ibyte + 1;   // (IBYTE + JBYTE + 1), since IBYTE==JBYTE for self
    int icont_off  = 1 + ipost;

    bzero(unit->buck, sizeof(int64) * NCONTS);

    uint8 *buf = malloc(1 << 20);
    if (buf == NULL) { perror("malloc recompute buf"); MPI_Abort(MPI_COMM_WORLD, 1); }
    int64 partial = 0;
    while (1) {
        ssize_t n = read(fd, buf + partial, (1 << 20) - partial);
        if (n < 0) { perror("read recompute"); MPI_Abort(MPI_COMM_WORLD, 1); }
        if (n == 0 && partial == 0) break;
        int64 avail = partial + n;
        int64 i = 0;
        while (i + seed_bytes <= avail) {
            int64 icont = 0;
            memcpy(&icont, buf + i + icont_off, icont_bytes);
            // strip sign bit (it lives in the high bit of the cont field)
            icont &= (((int64) 1 << (8 * icont_bytes - 1)) - 1);
            unit->buck[icont]++;
            i += seed_bytes;
        }
        partial = avail - i;
        if (partial > 0) memmove(buf, buf + i, partial);
        if (n == 0) break;
    }
    free(buf);
    (void) size;

    // Rewind so reimport_thread reads from start.
    if (lseek(fd, 0, SEEK_SET) < 0) { perror("lseek rewind"); MPI_Abort(MPI_COMM_WORLD, 1); }
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
    t_stage_start = MPI_Wtime();
    RecordBuffer local;
    record_buffer_init(&local, sizing);
    extract_to_records(gdb, world_size, my_rank, &local);
    t_extract = MPI_Wtime() - t_stage_start;

    // ---- STAGE: sort ----------------------------------------------------
    t_stage_start = MPI_Wtime();
    RecordBuffer recv;
    memset(&recv, 0, sizeof(recv));
    run_stage2(&local, &recv, MPI_COMM_WORLD);
    record_buffer_free(&local);
    t_sort = MPI_Wtime() - t_stage_start;

    fprintf(stderr, "[merge] rank %d: own slice has %lld records (sorted)\n",
            my_rank, (long long) recv.count);

    // ---- STAGE: merge ---------------------------------------------------
    // Covers: build_kmer_adapter + extern setup + IOBuffer alloc +
    // self_adaptamer_merge + redistribute_seeds(N+C) + recompute_buck.
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

    // Set FastGA externs (mirrors pga-merge.c).
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

    int64 total_buf_bytes = 2ll * num_units * 1000000ll;
    uint8 *unit_buffer = malloc(total_buf_bytes);
    int64 *unit_bucks  = calloc(2ll * NTHREADS * NCONTS, sizeof(int64));
    if (unit_buffer == NULL || unit_bucks == NULL) {
        fprintf(stderr, "rank %d: out of memory for unit buffers\n", my_rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int k = 0;
    for (int i = 0; i < NTHREADS; i++) {
        for (int j = 0; j < NUNITS; j++) {
            char npath[1024], cpath[1024];
            snprintf(npath, sizeof(npath), "%s/%s.%d.N", SORT_PATH, PAIR_NAME, k);
            snprintf(cpath, sizeof(cpath), "%s/%s.%d.C", SORT_PATH, PAIR_NAME, k);
            N_Units[k].bufr = unit_buffer + (2ll * k) * 1000000ll;
            C_Units[k].bufr = unit_buffer + (2ll * k + 1) * 1000000ll;
            N_Units[k].btop = N_Units[k].bufr;
            C_Units[k].btop = C_Units[k].bufr;
            N_Units[k].bend = N_Units[k].bufr + (1000000 - (IBYTE + JBYTE + 1));
            C_Units[k].bend = C_Units[k].bufr + (1000000 - (IBYTE + JBYTE + 1));
            N_Units[k].buck = unit_bucks + (2ll * i) * NCONTS;
            C_Units[k].buck = unit_bucks + (2ll * i + 1) * NCONTS;
            N_Units[k].inum = k;
            C_Units[k].inum = k;
            N_Units[k].file = open(npath, O_RDWR | O_CREAT | O_TRUNC, 0666);
            C_Units[k].file = open(cpath, O_RDWR | O_CREAT | O_TRUNC, 0666);
            if (N_Units[k].file < 0 || C_Units[k].file < 0) {
                fprintf(stderr, "rank %d: cannot open seed files\n", my_rank);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            k++;
        }
    }

    if (my_rank == ROOT)
        fprintf(stderr, "[merge] starting self-adaptamer merge on each rank...\n");

    // self_adaptamer_merge expects gdb->seqtot but we want the per-rank
    // record count proxy here for the verbose output ratio. Pass seqtot
    // since it's only used for a logging division.
    self_adaptamer_merge((Kmer_Stream *)&T1, (Post_List *)&P1, gdb->seqtot);

    // After the two Alltoallvs each rank's N_Units[0].file / C_Units[0].file
    // holds all seed pairs destined for it (from all source ranks).
    if (my_rank == ROOT)
        fprintf(stderr, "[merge] redistributing seeds via MPI_Alltoallv...\n");

    char n_combined[1024], c_combined[1024];
    int64 n_recv = redistribute_seeds(N_Units, world_size, 'N',
                                      SORT_PATH, PAIR_NAME,
                                      n_combined, sizeof(n_combined));
    int64 c_recv = redistribute_seeds(C_Units, world_size, 'C',
                                      SORT_PATH, PAIR_NAME,
                                      c_combined, sizeof(c_combined));
    fprintf(stderr, "[merge] rank %d: received N=%lld C=%lld bytes\n",
            my_rank, (long long) n_recv, (long long) c_recv);

    // After exchange we have ONE slot per type. NUNITS for downstream code
    // is now 1 (matching pair_sort_search's NPARTS=1 expectation).
    NUNITS = 1;

    // Recompute per-A-contig seed counts in the received buffers.
    recompute_buck(&N_Units[0], IBYTE, IPOST, ICONT);
    recompute_buck(&C_Units[0], IBYTE, IPOST, ICONT);
    t_merge = MPI_Wtime() - t_stage_start;

    // ---- STAGE: alignments ----------------------------------------------
    // Covers pair_sort_search + ALNtoPAF.
    t_stage_start = MPI_Wtime();

    // pair_sort_search uses NPARTS=1, so it only looks at N_Units[0] and
    // C_Units[0] — exactly the slots we filled with received bytes.
    fprintf(stderr, "[merge] rank %d: starting pair_sort_search\n", my_rank);
    pair_sort_search(gdb, gdb);

    char aln_path[1024];
    snprintf(aln_path, sizeof(aln_path), "%s/%s.1aln", SORT_PATH, ALGN_UNIQ);
    char cmd[2048];
    snprintf(cmd, sizeof(cmd),
             "ALNtoPAF %s > out.rank%d.paf 2>/dev/null", aln_path, my_rank);
    fprintf(stderr, "[merge] rank %d: %s\n", my_rank, cmd);
    int ret = system(cmd);
    if (ret != 0) fprintf(stderr, "rank %d: ALNtoPAF failed (ret=%d)\n", my_rank, ret);
    t_alignments = MPI_Wtime() - t_stage_start;

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

    char rm_cmd[1024];
    snprintf(rm_cmd, sizeof(rm_cmd), "rm -rf %s", SORT_PATH);
    system(rm_cmd);

    free(SORT_PATH);
    free(PAIR_NAME);
    free(ALGN_UNIQ);
    free(ALGN_PAIR);
    free(Perm1);
    free(IDBsplit);
    free(Select);
    free(unit_buffer);
    free(unit_bucks);
    free(N_Units);
    free(C_Units);

    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == ROOT) fprintf(stderr, "[merge] done\n");
    MPI_Finalize();
    return 0;
}
