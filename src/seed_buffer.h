#ifndef SEED_BUFFER_H
#define SEED_BUFFER_H

#include "gene_core.h"

typedef struct {
    uint8  *data;
    int64   len;
    int64   cap;
    int     inum;
} SeedUnit;

void seed_unit_init(SeedUnit *u, int inum);
void seed_unit_grow(SeedUnit *u, int64 need);
void seed_unit_trim(SeedUnit *u);
void seed_unit_free(SeedUnit *u);

static inline void seed_unit_ensure(SeedUnit *u, int64 need)
{ if (u->len + need > u->cap)
    seed_unit_grow(u, need);
}

#endif
