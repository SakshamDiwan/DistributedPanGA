#ifndef DPGA_SORT_H
#define DPGA_SORT_H

#include "extract.h"

// Sort a RecordBuffer in-place using msd_sort.
// Sorts by k-mer key (bytes 1..KBYTES of each record).
// After return, buf->data is in lex-ascending order on the k-mer key,
// and the LCP byte (byte 0) of each record holds msd_sort's LCP value.
void sort_records(RecordBuffer *buf);

#endif
