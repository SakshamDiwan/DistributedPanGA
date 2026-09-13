#ifndef PGA_SEED_EXCHANGE_H
#define PGA_SEED_EXCHANGE_H

#include <stdint.h>

// Plan an MPI_Alltoallv byte exchange: turn per-destination byte sizes into the
// int counts[] and displs[] arrays MPI requires.
//
// MPI-3 Alltoallv takes int counts and displacements, so both the individual
// sizes and the running displacement must fit in an int.
//
// Returns 0 on success. Negative on failure:
//   -1  a negative size
//   -2  a single destination exceeds INT_MAX
//   -3  the cumulative total exceeds INT_MAX (a later displacement would wrap)
//
// MPI-free by design, so the caller keeps the MPI_Abort and this stays
// trivially testable with large sizes and no allocation.
int plan_byte_exchange(const int64_t *sizes, int n, int *counts, int *displs);

#endif
