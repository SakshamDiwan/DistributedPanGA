// pga-prefix-hist: offline k-mer prefix histogram (R3 analysis tool).
//
// Runs the production extraction scan (extract_scan, the same k-mers as
// extract_to_records) over a GDB and counts, without storing any record:
//   <out>.hist24.bin   int64 little-endian [2^24]: records per 12-base prefix
//                      (first 3 k-mer bytes). The 1024 dsort buckets are the
//                      top 10 bits of this index (bucket = prefix24 >> 14).
//   <out>.contigs.tsv  contig, clen, records   (one line per scanned contig)
//
// Usage: pga-prefix-hist <gdb-stem> <out-prefix> [nchunks chunk]
// With nchunks/chunk, scans only work-split slice `chunk` of `nchunks`
// (compute_work_split, i.e. what rank `chunk` of `nchunks` extracts), so the
// whole GDB can be covered by several processes and the outputs summed.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "extract.h"
#include "constants.h"
#include "GDB.h"

typedef struct {
    int64_t *hist24;
    int64_t *ctg_records;
} Counts;

static void count_cb(void *user, int contig, int64_t position, int strand,
                     const uint8_t *kmer)
{
    Counts *c = (Counts *) user;
    (void) position; (void) strand;
    c->hist24[((uint32_t) kmer[0] << 16) | ((uint32_t) kmer[1] << 8) | kmer[2]] += 1;
    c->ctg_records[contig] += 1;
}

int main(int argc, char *argv[])
{
    if (argc != 3 && argc != 5) {
        fprintf(stderr, "Usage: %s <gdb-stem> <out-prefix> [nchunks chunk]\n", argv[0]);
        return 1;
    }
    int nchunks = 1, chunk = 0;
    if (argc == 5) {
        nchunks = atoi(argv[3]);
        chunk   = atoi(argv[4]);
        if (nchunks < 1 || chunk < 0 || chunk >= nchunks) {
            fprintf(stderr, "bad chunk %d of %d\n", chunk, nchunks);
            return 1;
        }
    }

    GDB _gdb, *gdb = &_gdb;
    // Read_GDB returns 1 (not < 0) on failure; check both.
    if (Read_GDB(gdb, argv[1]) != 0) {
        fprintf(stderr, "Read_GDB failed for '%s'\n", argv[1]);
        return 1;
    }

    Counts c;
    c.hist24      = (int64_t *) calloc((size_t) 1 << 24, sizeof(int64_t));
    c.ctg_records = (int64_t *) calloc(gdb->ncontig, sizeof(int64_t));
    if (c.hist24 == NULL || c.ctg_records == NULL) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }

    extract_stats_reset();
    extract_scan(gdb, nchunks, chunk, count_cb, &c);

    char path[4096];
    snprintf(path, sizeof(path), "%s.hist24.bin", argv[2]);
    FILE *f = fopen(path, "wb");
    if (f == NULL || fwrite(c.hist24, sizeof(int64_t), (size_t) 1 << 24, f) != ((size_t) 1 << 24)) {
        fprintf(stderr, "cannot write %s\n", path);
        return 1;
    }
    fclose(f);

    snprintf(path, sizeof(path), "%s.contigs.tsv", argv[2]);
    f = fopen(path, "w");
    if (f == NULL) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    int64_t total = 0;
    for (int r = 0; r < gdb->ncontig; r++) {
        total += c.ctg_records[r];
        if (c.ctg_records[r] > 0 || gdb->contigs[r].clen > 0)
            fprintf(f, "%d\t%lld\t%lld\n", r, (long long) gdb->contigs[r].clen,
                    (long long) c.ctg_records[r]);
    }
    fclose(f);

    fprintf(stderr, "pga-prefix-hist: %s chunk %d/%d: contigs=%lld bases=%lld records=%lld\n",
            argv[1], chunk, nchunks, (long long) g_extract_stats.contigs,
            (long long) g_extract_stats.bases, (long long) total);
    Close_GDB(gdb);
    return 0;
}
