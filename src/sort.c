#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "sort.h"
#include "record.h"
#include "constants.h"
#include "gene_core.h"          // brings in int64/uint8 used by msd_sort

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

    int64 bucket_start[256];
    int64 cum = 0;
    for (int i = 0; i < 256; i++)
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
