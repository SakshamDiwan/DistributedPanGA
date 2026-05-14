#ifndef CONTIG_ASSIGNMENT_H
#define CONTIG_ASSIGNMENT_H

#include <stdint.h>
#include "GDB.h"

// Deterministic contig-pair → rank assignment for self-alignment.
//
// For N contigs we have N*(N+1)/2 unordered pairs (i, j) with i ≤ j.
// Each pair gets assigned to exactly one rank in [0, world_size).
//
// The assignment is computed redundantly on every rank from GDB metadata
// (which all ranks already have), so no communication is needed.
//
// Policy (v1): weighted round-robin. Sort pairs by clen[i]*clen[j] descending,
// assign each to the currently least-loaded rank by accumulated weight.
// Stable / deterministic given identical (gdb, world_size) on every rank.

typedef struct {
    int      world_size;
    int      ncontig;
    // For pair (i, j) with i ≤ j, owner = rank_for_pair[i * ncontig + j].
    // Entries with i > j are unused (left as -1).
    int8_t  *rank_for_pair;
    // Per-rank load (sum of clen[i]*clen[j] for owned pairs). For diagnostics.
    int64_t *load_per_rank;
} ContigAssignment;

// Returns 0 on success, nonzero on alloc failure.
int  contig_assignment_build(const GDB *gdb, int world_size, ContigAssignment *out);

void contig_assignment_free(ContigAssignment *a);

static inline int contig_assignment_owner(const ContigAssignment *a, int i, int j)
{
    int lo = (i < j) ? i : j;
    int hi = (i < j) ? j : i;
    return a->rank_for_pair[(int64_t) lo * a->ncontig + hi];
}

void contig_assignment_print(const ContigAssignment *a, const GDB *gdb);

#endif
