#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "sort.h"
#include "record.h"
#include "constants.h"
#include "gene_core.h"          // brings in int64/uint8 used by msd_sort
#include <time.h>

static double sort_now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

extern void msd_sort(uint8 *array, int64 nelem, int rsize, int ksize,
                     int64 *part, int beg, int end, int nthreads);

// "Bases in common" between two k-mer byte arrays (2-bit packed).
// Returns count in [0, kbytes*4]. Used by recalc_all_lcps.
int simple_kmer_lcp_bases(const uint8_t *a, const uint8_t *b, int kbytes)
{
    int total_bases = kbytes * 4;
    int lcp = 0;
    for (int byte_idx = 0; byte_idx < kbytes; byte_idx++) {
        if (a[byte_idx] != b[byte_idx]) {
            uint8_t diff = a[byte_idx] ^ b[byte_idx];
            for (int base_in_byte = 0; base_in_byte < 4; base_in_byte++) {
                int shift = (3 - base_in_byte) * 2;
                if ((diff >> shift) & 0x03) return lcp;
                lcp++;
            }
        }
        lcp += 4;
    }
    return total_bases;
}

// Walk a sorted RecordBuffer and overwrite byte 0 of every record with the
// "bases in common" LCP relative to its predecessor. First record gets 0.
// Required because msd_sort's internal byte 0 encoding diverges from
// "bases in common" at part[] boundaries — see dsort.c §LCP encoding.
void recalc_all_lcps(RecordBuffer *buf)
{
    // An empty buffer is a genuine no-op, but a single record is not: it has no
    // predecessor, so its LCP is 0 -- exactly as for the first record of a
    // longer buffer. Returning early for count == 1 left whatever msd_sort had
    // written in byte 0. Reachable whenever a rank receives a one-record slice.
    if (buf->count == 0) return;
    int rsize = buf->sizing.record_size;
    buf->data[0] = 0;
    if (buf->count == 1) return;
    for (int64_t i = 1; i < buf->count; i++) {
        const uint8_t *prev = buf->data + (i - 1) * rsize + KMER_OFFSET;
        const uint8_t *curr = buf->data + i * rsize + KMER_OFFSET;
#ifdef LCPs
        buf->data[i * rsize] = (uint8_t) simple_kmer_lcp_bases(prev, curr, KBYTES);
#else
        buf->data[i * rsize] = (memcmp(prev, curr, KBYTES) == 0) ? 0 : 1;
#endif
    }
}

SortTiming g_sort_timing;

// Thread count and part granularity for sort_records. Defaults: 4 threads,
// 1024 parts (first 5 bases = the dsort bucket = the parts FastGA's GIXmake
// hands msd_sort). The environment overrides them:
//   PGA_SORT_THREADS = 1..256 (best: the rank's physical cores)
//   PGA_SORT_PARTS   = 256 (first k-mer byte; the pre-R3 behaviour) or 1024.
// Only part boundaries change: each part still holds whole runs of equal
// leading bytes, so the concatenation of the sorted parts is sorted.
int sort_threads_from_env(void)
{
    const char *e = getenv("PGA_SORT_THREADS");
    int v = (e != NULL && e[0] != '\0') ? atoi(e) : 0;
    return (v >= 1 && v <= 256) ? v : 4;
}

int sort_parts_from_env(void)
{
    const char *e = getenv("PGA_SORT_PARTS");
    int v = (e != NULL && e[0] != '\0') ? atoi(e) : 0;
    return (v == 256) ? 256 : 1024;
}

void sort_records(RecordBuffer *buf)
{
    sort_records_opts(buf, sort_threads_from_env(), sort_parts_from_env());
}

static inline int part_of(const uint8_t *rec, int nparts)
{
    // rec[1] is the first k-mer byte. 1024 parts = first 10 bits = the dsort
    // bucket (record_bucket in dsort.h).
    return (nparts == 1024) ? (((int) rec[1] << 2) | ((int) rec[2] >> 6)) : rec[1];
}

void sort_records_opts(RecordBuffer *buf, int nthreads, int nparts)
{
    int rsize = buf->sizing.record_size;
    // ksize = 1 (LCP byte) + KBYTES (k-mer bytes). msd_sort starts at digit=1.
    int ksize = 1 + KBYTES;
    double t0 = sort_now();

    if (nparts != 1024) nparts = 256;
    if (nthreads < 1) nthreads = 1;

    // Group the records by their leading bits (first five bases by default)
    // so we hand msd_sort a contiguous-by-part layout and it can sort the parts
    // in parallel, one thread per part. 1024 parts is FastGA's GIXmake layout.
    int64 bucket_count[1024];
    for (int i = 0; i < nparts; i++) bucket_count[i] = 0;
    for (int64 i = 0; i < buf->count; i++)
        bucket_count[part_of(buf->data + i * rsize, nparts)]++;

    int64 bucket_start[1024];
    int64 cum = 0;
    for (int i = 0; i < nparts; i++)
    {
        bucket_start[i] = cum * rsize;
        cum += bucket_count[i];
    }

    int64 total_bytes = (int64) buf->count * rsize;
    uint8_t *tmp = (uint8_t *) malloc(total_bytes + 1);
    if (tmp == NULL)
    {
        fprintf(stderr, "sort_records: out of memory for rebucket buffer\n");
        exit(1);
    }
    int64 cursor[1024];
    for (int i = 0; i < nparts; i++) cursor[i] = bucket_start[i];

    for (int64 i = 0; i < buf->count; i++)
    {
        uint8_t *src = buf->data + i * rsize;
        int      k0  = part_of(src, nparts);
        memcpy(tmp + cursor[k0], src, rsize);
        cursor[k0] += rsize;
    }
    memcpy(buf->data, tmp, total_bytes);
    free(tmp);

    // Build msd_sort's part[] array: byte length per part.
    // msd_sort writes a sentinel byte at array[total_bytes]; record_buffer
    // allocations include +1 for that.
    int64 part[1024];
    int64 largest = 0;
    int   nonempty = 0;
    for (int i = 0; i < nparts; i++)
    {
        part[i] = bucket_count[i] * rsize;
        if (bucket_count[i] > 0) nonempty++;
        if (bucket_count[i] > largest) largest = bucket_count[i];
    }
    double t1 = sort_now();

    // With 1024 parts msd_sort's own part-seam byte-0 values index
    // LCP_Table[x ^ n] only when parts x and n share bits 8-9, so x ^ n < 256
    // stays in range. Those byte-0 values are overwritten by recalc_all_lcps.
    msd_sort(buf->data, buf->count, rsize, ksize, part, 0, nparts, nthreads);

    g_sort_timing.rebucket_s     = t1 - t0;
    g_sort_timing.msd_s          = sort_now() - t1;
    g_sort_timing.threads        = nthreads;
    g_sort_timing.parts          = nparts;
    g_sort_timing.nonempty_parts = nonempty;
    g_sort_timing.largest_part   = largest;
}
