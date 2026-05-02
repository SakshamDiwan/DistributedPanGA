#ifndef DPGA_PACK_H
#define DPGA_PACK_H

#include <stdint.h>
#include "GDB.h"

// Pack the forward 40-mer at per-contig position j into 10 bytes.
// 40-mer covers forward bases [j, j+39].
// Returns 0 on success, -1 if [j, j+40) extends past contig.
int pack_forward_kmer(GDB *gdb, int contig_id, int64_t j, uint8_t *out);

// Pack the reverse-complement 40-mer keyed at scan-thread position j.
// Forward source covers bases [j-28, j+11] (reverse-complemented).
// Returns 0 on success, -1 if the source range extends past contig.
int pack_rc_kmer(GDB *gdb, int contig_id, int64_t j, uint8_t *out);

#endif
