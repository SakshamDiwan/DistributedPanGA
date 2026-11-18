#ifndef FASTGA_ORDER_H
#define FASTGA_ORDER_H

// Start each same-strand self pair from the value FastGA's alignment loop
// would have had on entry to it, so the output does not depend on the rank
// count. See fastga_order.c and Self_Order in fastga_pipeline.h.

#include <mpi.h>
#include <stdint.h>
#include "GDB.h"
#include "fastga_pipeline.h"

// This rank's same-strand seeds as received from the exchange: records of
// 1 + apos + actg + bpos + bctg bytes (length, A position, A contig, B
// position, B contig with the strand in its high bit).
typedef struct {
    const uint8_t *data;
    int64_t        len;
    int            apos_bytes, actg_bytes, bpos_bytes, bctg_bytes;
} FastgaOrderSeeds;

// Collective. Reads FastGA's contig order from <gdb_stem>.gix on rank 0 and
// installs Self_Order on every rank. n_seed_count[c] is this rank's number of
// same-strand seeds with anchor contig c; swide is the sort-array record
// width; rank_for_pair[a*ncontig+b] is the rank that owns contig pair (a,b).
// Returns 1 if enabled, 0 if the index could not be read (then
// nothing changes).
//   PGA_FASTGA_THREADS      thread count the reference FastGA run used (default 32)
//   PGA_FASTGA_UNSET_BEPOS  value where FastGA's main thread reads a field it
//                           never set (default 32767; other threads: 0)
//   PGA_FASTGA_ORDER_LOG    1: log every self pair's entry value and the ranges
int  fastga_order_setup(const char *gdb_stem, GDB *gdb, const int64_t *n_seed_count,
                        int swide, const FastgaOrderSeeds *seeds,
                        const int8_t *rank_for_pair, MPI_Comm comm);

// One summary line for this rank; uninstalls Self_Order.
void fastga_order_report(void);

#endif
