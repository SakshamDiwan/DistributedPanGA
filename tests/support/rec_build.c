#include <stdlib.h>
#include <string.h>
#include "rec_build.h"
#include "constants.h"

RecordSizing rs_make(int post_bytes, int cont_bytes)
{
    RecordSizing s;
    s.post_bytes  = post_bytes;
    s.cont_bytes  = cont_bytes;
    s.record_size = 1 + KBYTES + 1 + post_bytes + cont_bytes;
    return s;
}

// Allocate count == capacity == n at exactly n*record_size + 1 bytes.
//
// The +1 is msd_sort's sentinel slot -- it writes one byte past the last
// record.  Allocating with NO spare capacity is deliberate: the sentinel then
// lands on that final byte and nowhere else, so ASan catches a sorter that
// writes further.  A buffer with slack would silently absorb the overrun.
// See S-05a and S-05b, which split this concern.
void rb_alloc_exact(RecordBuffer *b, RecordSizing s, int64_t n)
{
    memset(b, 0, sizeof(*b));
    b->sizing   = s;
    b->count    = n;
    b->capacity = n;
    b->data     = (uint8_t *) calloc((size_t) (n * s.record_size + 1), 1);
}

void rb_free(RecordBuffer *b)
{
    free(b->data);
    memset(b, 0, sizeof(*b));
}

// Encode one record.  Note the two different byte orders in a single record:
// the k-mer is stored most-significant-base-first so memcmp sorts k-mers
// correctly, while position and contig are little-endian.  The strand flag is
// OR-ed into the top bit of the contig field.
void rb_set(RecordBuffer *b, int64_t i, const uint8_t *kmer,
            int64_t pos, int contig, int strand, uint8_t lcp)
{
    int      pb = b->sizing.post_bytes, cb = b->sizing.cont_bytes;
    uint8_t *r  = b->data + i * b->sizing.record_size;
    int64_t  cv = contig;
    int      k;

    r[LCP_OFFSET] = lcp;
    memcpy(r + KMER_OFFSET, kmer, KBYTES);
    r[MASK_OFFSET] = 0;
    for (k = 0; k < pb; k++) r[POS_OFFSET + k] = (uint8_t) ((pos >> (8 * k)) & 0xff);
    if (strand) cv |= ((int64_t) 1) << (8 * cb - 1);
    for (k = 0; k < cb; k++) r[POS_OFFSET + pb + k] = (uint8_t) ((cv >> (8 * k)) & 0xff);
}

// Decode one record back into plain fields, per the layout in rec_build.h.
// The strand flag occupies the top bit of the contig field, so the id is
// whatever remains below it -- hence masking rather than shifting.
void rb_decode(const RecordBuffer *b, int64_t i, RefRec *out, uint8_t *lcp_out)
{
    int            pb = b->sizing.post_bytes, cb = b->sizing.cont_bytes;
    const uint8_t *r  = b->data + i * b->sizing.record_size;
    int64_t        pos = 0, cv = 0, mask;
    int            k;

    for (k = 0; k < pb; k++) pos |= ((int64_t) r[POS_OFFSET + k])      << (8 * k);
    for (k = 0; k < cb; k++) cv  |= ((int64_t) r[POS_OFFSET + pb + k]) << (8 * k);

    mask = (((int64_t) 1) << (8 * cb - 1));
    out->position = pos;
    out->strand   = (cv & mask) ? 1 : 0;
    out->contig   = (int) (cv & (mask - 1));
    memcpy(out->kmer, r + KMER_OFFSET, KBYTES);
    if (lcp_out) *lcp_out = r[LCP_OFFSET];
}

void rb_decode_all(const RecordBuffer *b, RefRecList *out)
{
    int64_t i;
    ref_list_init(out);
    out->cap = b->count > 0 ? b->count : 1;
    out->v   = (RefRec *) malloc((size_t) out->cap * sizeof(RefRec));
    out->n   = b->count;
    for (i = 0; i < b->count; i++) rb_decode(b, i, &out->v[i], NULL);
}
