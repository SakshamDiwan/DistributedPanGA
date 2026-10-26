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

// --- Instrumentation (PR 1) -------------------------------------------------
// Phase times of the last pair_sort_search call, in seconds; index 0 = normal
// (N) seeds, 1 = complemented (C) seeds. Written by thread 0 only.
typedef struct {
    double reimport_s[2];   // reimport_thread: seed file -> sort array
    double sort_s[2];       // rmsd_sort
    double search_s[2];     // search_seeds (all align_contigs calls)
    double la_sort_s;
    double la_merge_s;
    int64  nels[2];         // seeds sorted per type
    int64  sarray_bytes;    // (nelmax+1)*swide, the sort array allocation
} PairSearchStats;

extern PairSearchStats Pair_Search_Stats;

// Optional per-task hook, called after every align_contigs call in
// search_seeds: one task = all seeds of one (anchor, partner) contig pair of
// one type. comp: 0 = N, 1 = C. nhits/nlass are the deltas for this task.
// NULL (the default) disables it. Must be thread-safe if NTHREADS > 1.
typedef void (*PairTaskHook)(int tid, int comp, int anchor, int partner,
                             int64 nseeds, double seconds,
                             int64 nhits, int64 nlass);
extern PairTaskHook Pair_Task_Hook;
// la_merge is internal to fastga_pipeline.c (called by pair_sort_search).

#endif