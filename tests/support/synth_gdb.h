#ifndef PGA_TEST_SYNTH_GDB_H
#define PGA_TEST_SYNTH_GDB_H

// Build a GDB entirely in memory from ACGT strings -- no files on disk, no
// FastGA tools, no gold fixtures.  This is what makes the extraction tests
// runnable anywhere in milliseconds.
//
// WHY THIS WORKS.  A GDB normally comes from Read_GDB reading a .1gdb plus a
// .bps file of 2-bit packed sequence.  But Get_Contig_Piece -- the only part
// extraction actually calls -- reads just four things from the struct:
// seqs, seqstate, ncontig, and contigs[i].{boff,clen}.  Everything else
// (provenance, scaffolds, headers, file paths) is untouched.  So we fill in
// those few fields by hand and point seqs at an in-memory stream.
//
// THE TRICK: seqstate = EXTERNAL makes Get_Contig_Piece treat seqs as a
// seekable FILE* over a 2-bit stream, so fmemopen() over a byte array works
// exactly like the real .bps file -- same fseeko, same fread, same
// Uncompress_Read.  This exercises byte-for-byte the code path production
// takes, since every driver reaches EXTERNAL via Read_GDB.
//
// DO NOT "simplify" THIS TO AN IN-MEMORY seqstate.  NUMERIC and COMPRESSED
// look easier -- seqs would just be a plain array -- but those branches of
// Get_Contig_Piece (lib/GDB.c:1868-1895) compute the source address as
// boff + beg/4 + beg%4, which is only correct for beg <= 3.  pack_forward_kmer
// is called with arbitrary offsets, so a NUMERIC fixture would silently return
// the wrong bases.  Production never reaches those branches, so the bug is
// harmless there -- but a test that used them would be chasing a phantom.
//
// OTHER RULES the fixture must respect:
//   - Get_Contig_Piece writes a sentinel at buffer[-1], so callers pass buf+1.
//   - Never call Close_GDB on one of these: it would free struct fields that
//     were never allocated.  synth_gdb_free releases only what we own.

#include <stdint.h>
#include <stdio.h>
#include "GDB.h"

typedef struct {
    GDB      gdb;
    uint8_t *bps;      // 2-bit packed blob backing gdb.seqs
    size_t   bps_len;
} SynthGDB;

// seqs[i] is a NUL-terminated ACGT string (upper case) for contig i.
// Returns 0 on success. Caller must synth_gdb_free().
int  synth_gdb_build(SynthGDB *sg, const char *const *seqs, int nseq);
void synth_gdb_free(SynthGDB *sg);

// ACGT -> 0..3, and back. 'A'=0 'C'=1 'G'=2 'T'=3.
int  base_code(char c);
char base_char(int v);

// Deterministic pseudorandom ACGT string of length len (xorshift, fixed seed).
// Caller frees.
char *synth_random_seq(int len, uint32_t seed);

#endif
