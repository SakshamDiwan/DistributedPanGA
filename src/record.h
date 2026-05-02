#ifndef DPGA_RECORD_H
#define DPGA_RECORD_H

#include <stdint.h>
#include "GDB.h"
#include "constants.h"

// Record layout (FastGA-compatible, except we store ALL 10 k-mer bytes
// instead of FastGA's 9 bytes with implicit-first-byte trick):
//
//   Byte 0:                          LCP placeholder (msd_sort fills this)
//   Bytes 1..KBYTES:                 packed 40-mer (10 bytes), MSB-first within byte
//   Byte KBYTES+1:                   mask placeholder (always 0; we don't support masking yet)
//   Bytes KBYTES+2..+1+post_bytes:   position within contig, little-endian
//   Bytes ...+1..+1+cont_bytes:      contig_id with strand bit in high position of last byte
//
// Total record size = 1 + KBYTES + 1 + post_bytes + cont_bytes.
// For yeast (longest contig ~1.5M, ~120 contigs): post_bytes=3, cont_bytes=1, record_size=16.

#define LCP_OFFSET     0
#define KMER_OFFSET    1
#define MASK_OFFSET    (1 + KBYTES)
#define POS_OFFSET     (1 + KBYTES + 1)

typedef struct {
    int post_bytes;
    int cont_bytes;
    int record_size;
} RecordSizing;

// Compute PostBytes/ContBytes from GDB metadata.
RecordSizing compute_record_sizing(GDB *gdb);

#endif
