#ifndef FASTGA_PIPELINE_H
#define FASTGA_PIPELINE_H

#include <stdint.h>
#include "GDB.h"
#include "libfastk.h"

// Post_List must match FastGA.c's definition (field order matters)
typedef struct {
    int     pbyte;
    int     cbyte;
    int64   nels;
    int64   maxp;
    int     freq;
    int     nctg;
    int    *perm;
    int64   cidx;
    uint8  *cache;
    uint8  *cptr;
    int64  *index;
    int     copn;
    int     part;
    int     nthr;
    int     nlen;
    char   *name;
    uint8  *ctop;
    int64  *neps;
    int     clone;
} Post_List;

// Per-thread / per-part IO buffer for seed pair files. Allocated by the
// driver (pga-merge.c), consumed by self_adaptamer_merge and pair_sort_search.
typedef struct {
    uint8  *bufr;
    uint8  *btop;
    uint8  *bend;
    int64  *buck;
    int     file;
    int     inum;
} IOBuffer;


void self_adaptamer_merge(Kmer_Stream *T1, Post_List *P1, int64 g1len);
void pair_sort_search(GDB *gdb1, GDB *gdb2);
// la_merge is internal to fastga_pipeline.c (called by pair_sort_search).

#endif