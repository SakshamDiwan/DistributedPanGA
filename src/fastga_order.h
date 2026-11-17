#ifndef FASTGA_ORDER_H
#define FASTGA_ORDER_H

// Start each same-strand self pair from the value FastGA's alignment loop
// would have had on entry to it, so the output does not depend on the rank
// count. See fastga_order.c and Self_Order in fastga_pipeline.h.

#include <mpi.h>
#include <stdint.h>
#include "GDB.h"
#include "fastga_pipeline.h"

// Collective. Reads FastGA's contig order from <gdb_stem>.gix on rank 0 and
// installs Self_Order on every rank. n_seed_count[c] is this rank's number of
// same-strand seeds with anchor contig c; swide is the sort-array record
// width. Returns 1 if enabled, 0 if the index could not be read (then
// nothing changes).
//   PGA_FASTGA_THREADS      thread count the reference FastGA run used (default 32)
//   PGA_FASTGA_UNSET_BEPOS  value where FastGA reads a field it never set (default 32767)
//   PGA_FASTGA_ORDER_LOG    1: log every self pair's entry value and the ranges
int  fastga_order_setup(const char *gdb_stem, GDB *gdb, const int64_t *n_seed_count,
                        int swide, MPI_Comm comm);

// One summary line for this rank; uninstalls Self_Order.
void fastga_order_report(void);

#endif
