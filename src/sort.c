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

int sort_heavy_from_env(void)
{
    const char *e = getenv("PGA_SORT_HEAVY");
    return !(e != NULL && strcmp(e, "0") == 0);
}

// Index of the next 8 key bits after the part prefix (10 bits for 1024
// parts, 8 for 256): a monotone function of the key within one part.
static inline int sub_of(const uint8_t *rec, int nparts)
{
    return (nparts == 1024) ? ((((int) rec[2] & 0x3f) << 2) | ((int) rec[3] >> 6)) : rec[2];
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

    // Heavy parts. msd_sort sorts each part on one thread, so a part holding
    // more than 1/nthreads of the records (poly-A/T buckets: up to 100 M
    // records on hsapiens90) bounds the sort. Each heavy part is regrouped
    // (stable) into 256 sub-parts on the next 8 key bits, and msd_sort gets
    // those sub-parts as separate parts. Part and sub-part boundaries are key
    // boundaries, so the result is sorted; only byte 0 (rewritten by
    // recalc_all_lcps) and the order of equal k-mers can differ.
    // PGA_SORT_HEAVY=0 disables this.
    int heavy_on = sort_heavy_from_env() && nthreads > 1;
    int64 heavy_min = buf->count / nthreads;
    if (heavy_min < 65536) heavy_min = 65536;
    int nheavy = 0;
    double t_heavy = 0;
    if (heavy_on)
        for (int i = 0; i < nparts; i++)
            if (bucket_count[i] > heavy_min) nheavy++;

    if (nheavy == 0)
    {
        // With 1024 parts msd_sort's own part-seam byte-0 values index
        // LCP_Table[x ^ n] only when parts x and n share bits 8-9, so x ^ n < 256
        // stays in range. Those byte-0 values are overwritten by recalc_all_lcps.
        msd_sort(buf->data, buf->count, rsize, ksize, part, 0, nparts, nthreads);
    }
    else
    {
        double th = sort_now();
        int64 maxheavy = 0;
        for (int i = 0; i < nparts; i++)
            if (bucket_count[i] > heavy_min && bucket_count[i] > maxheavy) maxheavy = bucket_count[i];
        int64 *heavy_sub[1024];
        for (int i = 0; i < nparts; i++) heavy_sub[i] = NULL;
        uint8_t *htmp = (uint8_t *) malloc(maxheavy * rsize + 1);
        if (htmp == NULL)
        {
            fprintf(stderr, "sort_records: out of memory for heavy-part buffer\n");
            exit(1);
        }
        for (int i = 0; i < nparts; i++)
        {
            if (bucket_count[i] <= heavy_min) continue;
            uint8_t *base = buf->data + bucket_start[i];
            int64 n = bucket_count[i];
            int64 sub_count[256], sub_cur[256];
            for (int j = 0; j < 256; j++) sub_count[j] = 0;
            for (int64 r = 0; r < n; r++)
                sub_count[sub_of(base + r * rsize, nparts)]++;
            int64 c = 0;
            for (int j = 0; j < 256; j++) { sub_cur[j] = c * rsize; c += sub_count[j]; }
            for (int64 r = 0; r < n; r++)
            {
                uint8_t *src = base + r * rsize;
                int j = sub_of(src, nparts);
                memcpy(htmp + sub_cur[j], src, rsize);
                sub_cur[j] += rsize;
            }
            memcpy(base, htmp, n * rsize);
            heavy_sub[i] = (int64 *) malloc(256 * sizeof(int64));
            if (heavy_sub[i] == NULL)
            {
                fprintf(stderr, "sort_records: out of memory for heavy sub-counts\n");
                exit(1);
            }
            memcpy(heavy_sub[i], sub_count, 256 * sizeof(int64));
        }
        free(htmp);
        // One msd_sort call over a combined list of contiguous key ranges: the
        // light parts and the heavy parts' sub-parts, in key order, merged
        // greedily to at most ~1000 segments (two consecutive segments always
        // exceed `target`, so there are at most 2*count/target + 1 = 1001).
        // A single call matters: msd_sort re-initialises its static mutex and
        // condition variable on entry while a worker of the previous call may
        // still be signalling, so back-to-back calls can lose a wakeup and hang.
        // Segment indices stay < 1024, so its seam code stays inside LCP_Table.
        int64 target = (buf->count + 499) / 500;     // ceil: count / target <= 500
        if (target < 1) target = 1;
        int64 *seg = (int64 *) malloc(1024 * sizeof(int64));
        if (seg == NULL)
        {
            fprintf(stderr, "sort_records: out of memory for segments\n");
            exit(1);
        }
        int nseg = 0;
        int64 cur = 0;                 // records in the open segment
#define ADD_UNIT(nrec) do { int64 u_ = (nrec);                                  \
            if (u_ == 0) break;                                               \
            if (cur > 0 && cur + u_ > target) { seg[nseg++] = cur * rsize; cur = 0; } \
            cur += u_; } while (0)
        for (int i = 0; i < nparts; i++)
        {
            if (bucket_count[i] <= heavy_min) { ADD_UNIT(bucket_count[i]); continue; }
            for (int j = 0; j < 256; j++) ADD_UNIT(heavy_sub[i][j]);
        }
        if (cur > 0) seg[nseg++] = cur * rsize;
#undef ADD_UNIT
        if (nseg > 1024)
        {
            fprintf(stderr, "sort_records: %d segments > 1024 (internal error)\n", nseg);
            exit(1);
        }
        msd_sort(buf->data, buf->count, rsize, ksize, seg, 0, nseg, nthreads);
        free(seg);
        for (int i = 0; i < nparts; i++) free(heavy_sub[i]);
        t_heavy = sort_now() - th;
    }

    g_sort_timing.rebucket_s     = t1 - t0;
    g_sort_timing.msd_s          = sort_now() - t1;
    g_sort_timing.heavy_parts    = nheavy;
    g_sort_timing.heavy_s        = t_heavy;
    g_sort_timing.threads        = nthreads;
    g_sort_timing.parts          = nparts;
    g_sort_timing.nonempty_parts = nonempty;
    g_sort_timing.largest_part   = largest;
}
