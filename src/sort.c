#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "sort.h"
#include "constants.h"
#include "gene_core.h"          // brings in int64/uint8 used by msd_sort

extern void msd_sort(uint8 *array, int64 nelem, int rsize, int ksize,
                     int64 *part, int beg, int end, int nthreads);

void sort_records(RecordBuffer *buf)
{
    int rsize = buf->sizing.record_size;
    // ksize = 1 (LCP byte) + KBYTES (k-mer bytes). msd_sort starts at digit=1.
    int ksize = 1 + KBYTES;

    // Bucket the records by first k-mer byte (byte 1) so we hand msd_sort
    // a contiguous-by-bucket layout and can let it parallelise per-bucket.
    // Even though for Phase 2 we'd be fine with one big partition, the
    // pre-bucketing matches FastGA's internal layout and lets the sort
    // multi-thread.
    int64 bucket_count[256];
    for (int i = 0; i < 256; i++) bucket_count[i] = 0;
    for (int64 i = 0; i < buf->count; i++)
    {
        uint8_t first = buf->data[i * rsize + 1];
        bucket_count[first]++;
    }

    // Prefix sum -> bucket starts (in bytes).
    int64 bucket_start[256];
    int64 cum = 0;
    for (int i = 0; i < 256; i++)
    {
        bucket_start[i] = cum * rsize;
        cum += bucket_count[i];
    }

    // Out-of-place rebucket (allocate a temp buffer).
    int64 total_bytes = (int64) buf->count * rsize;
    uint8_t *tmp = (uint8_t *) malloc(total_bytes + 1);
    if (tmp == NULL)
    {
        fprintf(stderr, "sort_records: out of memory for rebucket buffer\n");
        exit(1);
    }
    int64 cursor[256];
    for (int i = 0; i < 256; i++) cursor[i] = bucket_start[i];

    for (int64 i = 0; i < buf->count; i++)
    {
        uint8_t *src = buf->data + i * rsize;
        uint8_t  k0  = src[1];
        memcpy(tmp + cursor[k0], src, rsize);
        cursor[k0] += rsize;
    }
    memcpy(buf->data, tmp, total_bytes);
    free(tmp);

    // Build msd_sort's part[] array: byte length per first-byte bucket.
    // msd_sort writes a sentinel byte at array[total_bytes]; record_buffer
    // allocations include +1 for that.
    int64 part[256];
    for (int i = 0; i < 256; i++)
        part[i] = bucket_count[i] * rsize;

    int nthreads = 4;
    msd_sort(buf->data, buf->count, rsize, ksize, part, 0, 256, nthreads);
}
