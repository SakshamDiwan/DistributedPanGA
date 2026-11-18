#ifndef FASTGA_ORDER_CORE_H
#define FASTGA_ORDER_CORE_H

// MPI-free part of fastga_order.c: FastGA's contig order and thread ranges,
// and the rules that give a same-strand self pair the end position FastGA's
// alignment loop would have had on entry to it. Kept separate so it can be
// unit tested in one process (tests/unit/test_fastga_order.c).
//
// Positions x, y are positions in FastGA's contig order (the .gix
// permutation). A key packs a partner position and an end position so that
// MPI_MAX picks the latest partner: (pos+1) << 32 | bepos.

#include <stdint.h>

static inline int64_t fo_make_key(int pos_b, int bepos)
{ return (((int64_t) (pos_b+1)) << 32) | (uint32_t) bepos; }
static inline int fo_key_pos(int64_t key)   { return (int) (key >> 32) - 1; }
static inline int fo_key_bepos(int64_t key) { return (int) (key & 0xffffffffll); }

// Self pair states in S[x]: none, not yet aligned, aligned without a real
// aligner call, or FO_S_CALLED + final bepos.
#define FO_S_NONE    0
#define FO_S_PENDING 1
#define FO_S_NOCALL  2
#define FO_S_CALLED  3

// Answers about earlier partners, by nearness k (0 = nearest).
#define FO_ANSWER_UNKNOWN (-1)
#define FO_ANSWER_NOCALL  (-2)

// The contig permutation stored in <stem>.gix (layout: FastGA.c
// Open_Post_List). Returns 1 and fills perm[0..n) if the file has n contigs
// and a valid permutation; otherwise 0 and a reason in err.
int fo_read_gix_perm(const char *stem, int n, int *perm, char *err, int errlen);

// FastGA's split of the contigs into at most T parts by length (FastGA.c main,
// "Compute GDB split into NTHREADS parts") and of each part into at most T
// thread ranges by seed bytes (RSDsort.c rmsd_sort). clen[x] and seeds[x] are
// by position; npost is the GDB's total bases (gdb->seqtot, which FastGA
// uses); swide is the sort-array record width. Fills range_start[x]
// (first position of x's range; x itself where x is in no range) and, if beg
// and end are not NULL, the ranges [beg[i], end[i]) in order (room for n).
// Returns the number of ranges; *nparts gets the number of parts.
int fo_compute_ranges(int n, const int64_t *clen, int64_t npost, const int64_t *seeds,
                      int64_t swide, int T, int *range_start, int *beg, int *end, int *nparts);

// has[xa*n+xb] != 0 for xb < xa: pair (xa,xb) has same-strand seeds.
// nearness: how many such partners lie strictly between xb and xa.
// count_earlier: how many such partners xa has. kth_earlier: the position of
// the k-th nearest (k = 0 nearest), or -1.
int fo_nearness(const uint8_t *has, int n, int xa, int xb);
int fo_count_earlier(const uint8_t *has, int n, int xa);
int fo_kth_earlier(const uint8_t *has, int n, int xa, int k);

// The value FastGA reads where it has no earlier aligner call in a range: a
// variable it never set. Observed 32764..32767 in the run's first range and 0
// in the others; unset_main is the value used for the first range.
int fo_unset_value(int range_start_x, int unset_main);

// Entry value of self pair x from its earlier partners' answers.
//   1  known: *bepos, and *k_used = nearness of the partner it came from, or
//      -1 if x is first in its range and no earlier partner made a call
//      (unset value)
//   0  an answer it needs has not arrived
//  -1  no earlier partner made a call and x is not first in its range: the
//      value comes from an earlier row (fo_from_rows)
int fo_from_answers(const int *answer, int nnear, int x, int range_start_x, int unset_main,
                    int *bepos, int *k_used);

// Entry value of self pair x from earlier rows of its range, given for each
// position y the key of the latest partner of y (other than y) that made a
// call, L[y], and the state of self pair (y,y), S[y]. In row y the order is:
// partners before y, then (y,y), then partners after y.
//   1  known: *bepos, *row = the row it came from and *from_pos = the partner
//      position (y itself for the self pair), or *row = -1 for the unset value
//   0  a self pair it depends on has not been aligned yet
int fo_from_rows(const int64_t *L, const int64_t *S, int x, int range_start_x, int unset_main,
                 int *bepos, int *row, int *from_pos);

#endif
