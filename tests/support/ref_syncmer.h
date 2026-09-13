#ifndef PGA_TEST_REF_SYNCMER_H
#define PGA_TEST_REF_SYNCMER_H

// An INDEPENDENT reference implementation of the extraction contract.
//
// This is the oracle the extraction tests compare against, and its value comes
// entirely from being written separately.  It is deliberately naive: a plain
// loop over every offset, recomputing all five hashes each time, with none of
// the rolling state, ring buffers or incremental minimum tracking that make
// src/extract.c fast.  If both implementations agree on thousands of offsets,
// it is very unlikely both are wrong in the same way.
//
// WHAT IT SHARES WITH PRODUCTION, and why that is acceptable: only TMap[], a
// fixed 256-entry scramble table.  That table is part of the ALGORITHM'S
// DEFINITION -- like a magic constant in a hash function -- not part of its
// logic, and reproducing it here would just be copying 256 numbers.  The
// reverse complement, by contrast, IS logic, so ref_rc4() recomputes it from
// first principles instead of reusing production's Comp[] lookup.
//
// THE ALGORITHM: CLOSED (12,8)-SYNCMERS.
//
// A syncmer scheme decides which k-mers to keep by looking at the smaller
// substrings inside them, so that the decision is consistent wherever the
// k-mer appears.  Here the 12-mer at offset j contains five 8-mers, starting
// at offsets j, j+1, j+2, j+3, j+4:
//
//     offset:  j   j+1  j+2  j+3  j+4
//     8-mers: [--------]
//                  [--------]
//                       [--------]
//                            [--------]
//                                 [--------]
//              \_____________ 12-mer ______/
//
// Each 8-mer gets a CANONICAL hash: the smaller of its forward hash and its
// reverse-complement hash, so a sequence and its reverse complement select the
// same positions.  The 12-mer is a CLOSED syncmer -- and is kept -- if and only
// if the smallest of those five hashes lies at either END of the window, at
// offset j or j+4.  Ties count as selected.
//
// Each selected offset j can then emit up to two records:
//     forward 40-mer, if j <= len-40 : stored at position j,    strand 0
//     RC 40-mer,      if j >= 28     : stored at position j+12, strand 1

#include <stdint.h>

typedef struct {
    int     contig;
    int64_t position;
    int     strand;        // 0 = forward, 1 = reverse complement
    uint8_t kmer[10];
} RefRec;

typedef struct {
    RefRec *v;
    int64_t n, cap;
    // Anti-vacuity counters, so a test can prove each selection branch was hit.
    int64_t n_right;       // index 4 is the strict unique minimum
    int64_t n_left;        // earliest argmin is index 0
    int64_t n_tie;         // index 4 ties the minimum without being earliest
    int64_t n_sel;         // selected 12-mers
} RefRecList;

void ref_list_init(RefRecList *l);
void ref_list_free(RefRecList *l);

// Is the 12-mer at offset j a closed syncmer? Used to build fixtures that
// require a specific offset to be selected.
int ref_is_selected(const char *seq, int len, int64_t j);

// Append the records extraction must produce for one contig.
void ref_extract_contig(const char *seq, int len, int contig, RefRecList *out);

// Independent packers, from base strings.
void ref_pack_forward(const char *seq, int64_t j, uint8_t *out);  // bases [j, j+40)
void ref_pack_rc(const char *seq, int64_t j, uint8_t *out);       // RC of [j-28, j+12)

// Independent base-level LCP: unpack both k-mers to base strings and count the
// common prefix. Returns [0, 40]. Used instead of production
// simple_kmer_lcp_bases so LCP checks are not self-referential.
int ref_lcp_bases(const uint8_t *a, const uint8_t *b);

// Sort + compare helpers for occurrence-multiset equality.
int  ref_rec_cmp(const void *x, const void *y);
void ref_list_sort(RefRecList *l);

#endif
