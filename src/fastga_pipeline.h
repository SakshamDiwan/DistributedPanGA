#ifndef FASTGA_PIPELINE_H
#define FASTGA_PIPELINE_H

#include <stdint.h>
#include "GDB.h"
#include "libfastk.h"
#include "seed_buffer.h"
#include "kmer_adapter.h"   // KmerStreamAdapter, PayloadLayout

// Per-thread / per-part per-contig seed counts (buck). Allocated by
// pga-mpi-merge.c. (FastGA's file-buffer fields were unused here: seeds live
// in SeedUnit buffers instead.)
typedef struct {
    int64  *buck;
} IOBuffer;


void self_adaptamer_merge(KmerStreamAdapter *kmer_table, PayloadLayout *payload_layout, int64 g1len);
void pair_sort_search(GDB *gdb1, GDB *gdb2);

// --- Instrumentation (PR 1) -------------------------------------------------
// Phase times of the last pair_sort_search call, in seconds; index 0 = normal
// (N) seeds, 1 = complemented (C) seeds. Written by thread 0 only.
typedef struct {
    double reimport_s[2];   // reimport_thread: seed buffer -> sort array
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
// one type. NULL (the default) disables it. Must be thread-safe if
// NTHREADS > 1.
typedef struct {
    int    tid;
    int    comp;        // 0 = N (normal) seeds, 1 = C (complemented)
    int    anchor;      // contig ids as passed to align_contigs
    int    partner;
    int64  alen;        // anchor / partner contig lengths
    int64  blen;
    // Known before alignment (from the task's seeds, sorted by diagonal):
    int64  nseeds;
    int64  ndiag;       // distinct diagonals among the seeds
    int64  maxdiag;     // most seeds on any single diagonal
    // Measured:
    double seconds;     // align_contigs wall time
    int64  nhits;       // deltas of the Contig_Bundle counters for this task
    int64  nlass;       // alignments found
    int64  nlive;       // non-redundant alignments
    int64  nlcov;       // bases covered by the non-redundant alignments
} PairTaskRecord;

typedef void (*PairTaskHook)(const PairTaskRecord *rec);
extern PairTaskHook Pair_Task_Hook;
// --- FastGA-order emulation ---------------------------------------------------
// FastGA's alignment loop skips the aligner for a same-strand self hit that
// straddles the main diagonal and then reads path->bepos, which still holds the
// end of the previous real aligner call. FastGA aligns pairs in a fixed order,
// so that value is fixed there; here it would depend on which pairs this rank
// owns. With Self_Order set, each self pair starts from the value FastGA would
// have had. That value is the final bepos of another pair, usually the pair
// with the nearest earlier partner in FastGA's contig order (the "candidate"),
// so the same-strand pairs are aligned candidates first, then the other pairs,
// and self pairs last, once the candidates' results have been exchanged. A self
// pair whose value is not yet known is set aside and aligned after all other
// pairs. The hooks hold the cross-rank bookkeeping (src/fastga_order.c). NULL
// (the default) keeps the single pass in sorted order. Requires NTHREADS == 1.
typedef struct {
    int  (*pair_class)(int a, int b);              // -2 candidate, -1 other partner before the
                                                   //   anchor, 0 self, 1 partner after it
    void (*record)(int a, int b, int calls, int bepos);  // after every same-strand pair
    void (*pass_done)(int pass);                   // 0 candidates, 1 other earlier, 2 later partners
    void (*poll)(void);                            // called often: lets the exchanges progress
    int  (*self_value)(int c, int *bepos);         // 1: entry value known now; 0: set the pair aside
    void (*self_begin)(int stage);                 // all other pairs done: finish exchange 0, then 1
    int  (*self_ready)(int c, int *bepos);         // set-aside pair: 1 if an exchange gave the value
    void (*round_begin)(void);                     // remaining pairs: exchange what is known
    int  (*resolve)(int c, int *bepos);            // 1: entry value known now; 0: wait a round
    int  (*round_end)(int64_t remaining);          // 1 while any rank still has pairs set aside
} SelfOrderHooks;
extern SelfOrderHooks *Self_Order;

// la_merge is internal to fastga_pipeline.c (called by pair_sort_search).

#endif