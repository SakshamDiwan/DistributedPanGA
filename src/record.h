#ifndef DPGA_RECORD_H
#define DPGA_RECORD_H

#include <stdint.h>
#include "GDB.h"
#include "constants.h"

// Record layout (FastGA-compatible, except we store ALL 10 k-mer bytes
// instead of FastGA's 9 bytes with implicit-first-byte trick):
//
//   Byte  LCP_OFFSET (0):                 LCP: bases in common with the previous
//                                         record's k-mer, [0,40]. Written by
//                                         recalc_all_lcps (run_stage2), not msd_sort.
//   Bytes KMER_OFFSET .. KBYTES (1..10):  packed 40-mer, 4 bases per byte, MSB-first
//   Byte  MASK_OFFSET (11):               mask, always 0 (masking not supported)
//   Bytes POS_OFFSET .. POS_OFFSET+position_bytes-1:
//                                         position within contig, little-endian
//   Bytes POS_OFFSET+position_bytes .. record_size-1:
//                                         contig id, little-endian; the strand is
//                                         the top bit of the last byte
//
// record_size = 1 + KBYTES + 1 + position_bytes + contig_strand_bytes.
// yeast7 (longest contig 1,566,853, 119 contigs): position_bytes=3,
// contig_strand_bytes=1, record_size=16 -> position = bytes 12..14,
// contig+strand = byte 15.

#define LCP_OFFSET     0
#define KMER_OFFSET    1
#define MASK_OFFSET    (1 + KBYTES)
#define POS_OFFSET     (1 + KBYTES + 1)

typedef struct {
    int position_bytes;
    int contig_strand_bytes;
    int record_size;
} RecordSizing;

// Compute position_bytes, contig_strand_bytes and record_size from GDB metadata
// (gdb->maxctg and gdb->ncontig).
RecordSizing compute_record_sizing(GDB *gdb);

#endif
