#ifndef DPGA_SORT_H
#define DPGA_SORT_H

#include "extract.h"

// Sort a RecordBuffer in-place using msd_sort.
// Sorts by k-mer key (bytes 1..KBYTES of each record).
// After return, buf->data is in lex-ascending order on the k-mer key,
// and the LCP byte (byte 0) of each record holds msd_sort's LCP value.
void sort_records(RecordBuffer *buf);

// Same, with explicit thread count and part granularity (256 = first k-mer
// byte, the default; 1024 = first five bases). sort_records reads them from
// PGA_SORT_THREADS / PGA_SORT_PARTS (defaults 4 / 256).
void sort_records_opts(RecordBuffer *buf, int nthreads, int nparts);
int  sort_threads_from_env(void);
int  sort_parts_from_env(void);

// Phase times of the most recent sort_records call (seconds, this process).
typedef struct {
    double  rebucket_s;      // part counting + out-of-place regroup + copy back
    double  msd_s;           // msd_sort (parallel over parts)
    int     threads, parts;  // settings used
    int     nonempty_parts;  // parts holding records (= available parallel tasks)
    int64_t largest_part;    // records in the largest part
} SortTiming;
extern SortTiming g_sort_timing;

// Overwrite byte 0 of every record with "bases in common" LCP vs predecessor.
void recalc_all_lcps(RecordBuffer *buf);

// "Bases in common" helper, exposed so dsort.c (cross-rank fixup) can use it.
int simple_kmer_lcp_bases(const uint8_t *a, const uint8_t *b, int kbytes);

#endif
