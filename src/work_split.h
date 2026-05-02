#ifndef DPGA_WORK_SPLIT_H
#define DPGA_WORK_SPLIT_H

#include <stdint.h>
#include "GDB.h"

// Divide contigs across workers so each gets roughly equal total bases.
// Outputs:
//   split[0..num_workers]  contig-index boundaries; worker w owns [split[w], split[w+1])
//   post[0..num_workers]   cumulative base count where worker w's region starts
//                          (post[0] = 0, post[num_workers] = gdb->seqtot)
// Caller must allocate split[num_workers+1] and post[num_workers+1].
void compute_work_split(GDB *gdb, int num_workers, int *split, int64_t *post);

#endif
