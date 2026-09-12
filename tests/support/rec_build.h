#ifndef PGA_TEST_REC_BUILD_H
#define PGA_TEST_REC_BUILD_H

// Hand-build and independently DECODE RecordBuffers.
//
// Two jobs.  Encoding lets sort/LCP tests construct a buffer with exactly the
// k-mers they want, without going through extraction at all.  Decoding lets
// extraction tests read production's output back into comparable values.
//
// The decoder is written from the layout spec in src/record.h, NOT by calling
// any production decode path -- the only such path is inside the pga-mpi
// driver's text dump, and using it would make extraction tests circular.
//
// RECORD LAYOUT.  Fixed prefix, then two variable-width fields whose widths
// come from compute_record_sizing:
//
//   byte 0                      LCP (bases shared with the previous record)
//   bytes 1..10                 the packed 40-mer, 4 bases per byte
//   byte 11                     mask byte, always 0 (masking unsupported)
//   POS_OFFSET .. +post_bytes   position within the contig, LITTLE-endian
//   .. +cont_bytes              contig id, little-endian, with the STRAND
//                               flag in the high bit of the field
//
// Note the two different byte orders in one record: the k-mer packs bases
// most-significant-first so that memcmp sorts k-mers correctly, while the
// numeric fields are little-endian.  Mixing them up is an easy mistake and is
// exactly what E-10 and E-07 check.

#include <stdint.h>
#include "extract.h"
#include "record.h"
#include "ref_syncmer.h"

RecordSizing rs_make(int post_bytes, int cont_bytes);

// Allocate with count == capacity == n and exactly n*record_size + 1 bytes,
// so msd_sort's sentinel write lands on the +1 byte and nowhere else.
void rb_alloc_exact(RecordBuffer *b, RecordSizing s, int64_t n);
void rb_free(RecordBuffer *b);

void rb_set(RecordBuffer *b, int64_t i, const uint8_t *kmer,
            int64_t pos, int contig, int strand, uint8_t lcp);

// Decode record i. lcp_out may be NULL.
void rb_decode(const RecordBuffer *b, int64_t i, RefRec *out, uint8_t *lcp_out);

// Decode every record into a freshly allocated RefRecList.
void rb_decode_all(const RecordBuffer *b, RefRecList *out);

#endif
