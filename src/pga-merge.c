#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "record.h"
#include "extract.h"
#include "sort.h"
#include "constants.h"
#include "GDB.h"
#include "kmer_adapter.h"
#include "libfastk.h"
#include "align.h"
#include "alncode.h"
#include "gene_core.h"
#include "fastga_pipeline.h"

// constants.h defines KMER as the literal 40, but FastGA's pipeline uses
// KMER as a runtime int (set from T1.kmer). Drop the macro so the extern
// declaration and the assignment below resolve to the variable.
#undef KMER

// Externs consumed by fastga_pipeline.c.
int    NTHREADS   = 1;
int    NPARTS     = 1;
int    SELF       = 1;
int    SOFT_MASK  = 0;
int    VERBOSE    = 1;
int    KMER       = 40;

int    FREQ       = 10;
int    CHAIN_BREAK = 2000;
int    CHAIN_MIN  = 170;
int    ALIGN_MIN  = 100;
double ALIGN_RATE = 0.3;
// Prog_Name lives in lib/gene_core.c; we set it in main().

int    IBYTE, JBYTE;
int    ICONT, JCONT;
int    IPOST, JPOST;
int    ISIGN, JSIGN;
int    KBYTE, CBYTE, LBYTE, PAYOFF;
int    ESHIFT;

int    NCONTS;
int64  AMXPOS, BMXPOS, MAXDAG;
int    DBYTE;

// Phase 4b: per-rank seed routing. NULL/1 means single-rank mode (use Select).
int8_t *RankForPair = NULL;
int     NUNITS      = 1;

int   *Select   = NULL;
int   *IDBsplit = NULL;
int   *Perm1    = NULL;
int   *Perm2    = NULL;
char  *SORT_PATH = NULL;
char  *PAIR_NAME = NULL;
char  *ALGN_UNIQ = NULL;
char  *ALGN_PAIR = NULL;

IOBuffer *N_Units = NULL;
IOBuffer *C_Units = NULL;

static char *make_temp_dir(void)
{
    char template[] = "/tmp/pga-merge-XXXXXX";
    char *dir = mkdtemp(template);
    if (dir == NULL) { perror("mkdtemp"); exit(1); }
    return strdup(dir);
}

int main(int argc, char *argv[])
{
    Prog_Name = "pga-merge";

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <gdb-stem>\n", argv[0]);
        return 1;
    }

    char *gdb_stem = argv[1];

    GDB _gdb, *gdb = &_gdb;
    if (Read_GDB(gdb, gdb_stem) < 0) {
        fprintf(stderr, "Read_GDB failed for '%s'\n", gdb_stem);
        return 1;
    }
    fprintf(stderr, "GDB: %d contigs, %lld total bases\n",
            gdb->ncontig, (long long) gdb->seqtot);

    RecordSizing sizing = compute_record_sizing(gdb);
    fprintf(stderr, "Record size: %d bytes (post=%d, cont=%d)\n",
            sizing.record_size, sizing.post_bytes, sizing.cont_bytes);

    RecordBuffer buf;
    record_buffer_init(&buf, sizing);
    extract_to_records(gdb, 1, 0, &buf);
    fprintf(stderr, "Extracted %lld records\n", (long long) buf.count);

    sort_records(&buf);

    recalc_all_lcps(&buf);
    fprintf(stderr, "Sorted and LCP-recalculated %lld records\n", (long long) buf.count);

    KmerStreamAdapter T1;
    PostListAdapter P1;
    memset(&T1, 0, sizeof(T1));
    memset(&P1, 0, sizeof(P1));

    if (build_kmer_adapter(&buf, &T1, &P1) < 0) {
        fprintf(stderr, "build_kmer_adapter failed\n");
        return 1;
    }
    fprintf(stderr, "Adapter built: %lld entries, pbyte=%d\n",
            (long long) T1.nels, T1.pbyte);

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

    // Our extract_to_records writes GDB-native contig ids into the cont
    // field (extract.c:270-276), so the perm that maps "id-as-stored" to
    // "GDB-native id" is the identity. Stock FastGA uses GIXmake's
    // length-sorted permutation here.
    Perm1 = malloc(gdb->ncontig * sizeof(int));
    Perm2 = Perm1;
    for (int i = 0; i < gdb->ncontig; i++) Perm1[i] = i;

    NPARTS = 1;
    IDBsplit = malloc((NPARTS + 1) * sizeof(int));
    IDBsplit[0] = 0;
    IDBsplit[1] = gdb->ncontig;

    Select = malloc(gdb->ncontig * sizeof(int));
    for (int i = 0; i < gdb->ncontig; i++) Select[i] = 0;

    SORT_PATH = make_temp_dir();
    fprintf(stderr, "Temp dir: %s\n", SORT_PATH);

    PAIR_NAME = Malloc(256, "pair_name");
    ALGN_UNIQ = Malloc(256, "algn_uniq");
    ALGN_PAIR = Malloc(256, "algn_pair");
    sprintf(PAIR_NAME, "_pair.%d", getpid());
    sprintf(ALGN_UNIQ, "_uniq.%d", getpid());
    sprintf(ALGN_PAIR, "_algn.%d", getpid());

    int num_units = NTHREADS * NPARTS;
    N_Units = calloc(num_units, sizeof(IOBuffer));
    C_Units = calloc(num_units, sizeof(IOBuffer));

    // One contiguous buffer for all N + C bufrs, matching FastGA.c:5106.
    // reimport_thread reads `2*NPARTS*1000000` bytes into N_Units[p].bufr
    // (fastga_pipeline.c:1103-1106), so the N/C bufrs must be adjacent in
    // memory to absorb that overflow read.
    int64 total_buf_bytes = 2ll * num_units * 1000000ll;
    uint8 *unit_buffer = malloc(total_buf_bytes);
    int64 *unit_bucks = calloc(2ll * NTHREADS * NCONTS, sizeof(int64));
    if (unit_buffer == NULL || unit_bucks == NULL) {
        fprintf(stderr, "out of memory for unit buffers\n");
        return 1;
    }

    int k = 0;
    for (int i = 0; i < NTHREADS; i++) {
        for (int j = 0; j < NPARTS; j++) {
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
                fprintf(stderr, "Cannot open seed-pair files in %s\n", SORT_PATH);
                return 1;
            }
            k++;
        }
    }

    fprintf(stderr, "Starting self-adaptamer merge...\n");
    self_adaptamer_merge((Kmer_Stream *)&T1, (Post_List *)&P1, gdb->seqtot);
    fprintf(stderr, "Merge complete\n");

    fprintf(stderr, "Starting pair sort and alignment search...\n");
    pair_sort_search(gdb, gdb);
    fprintf(stderr, "Alignments complete\n");

    char aln_path[1024];
    snprintf(aln_path, sizeof(aln_path), "%s/%s.1aln", SORT_PATH, ALGN_UNIQ);

    char cmd[2048];
    snprintf(cmd, sizeof(cmd), "ALNtoPAF %s > out.paf 2>/dev/null", aln_path);
    fprintf(stderr, "Running: %s\n", cmd);
    int ret = system(cmd);
    if (ret != 0) fprintf(stderr, "ALNtoPAF failed (ret=%d)\n", ret);

    free_kmer_adapter(&T1, &P1);
    record_buffer_free(&buf);
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

    fprintf(stderr, "Done. Output: out.paf\n");
    return 0;
}