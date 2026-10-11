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

// The same packings from bases already in memory (NUMERIC, one base per byte).
// Forward: bases[0..39] are contig bases [j, j+40).
// RC:      bases[0..39] are contig bases [j-28, j+12) (forward strand).
// pack_*_kmer above fetch the bases with Get_Contig_Piece and then call these.
void pack_forward_bases(const uint8_t *bases, uint8_t *out);
void pack_rc_bases(const uint8_t *bases, uint8_t *out);

#endif
