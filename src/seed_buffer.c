#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "seed_buffer.h"

#define SEED_UNIT_INIT_CAP  (1024 * 1024)

void seed_unit_init(SeedUnit *u, int inum)
{ u->data = malloc(SEED_UNIT_INIT_CAP);
  if (u->data == NULL)
    { int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      fprintf(stderr, "[rank %d] seed_unit_init: malloc(%d) failed for dest %d\n",
              rank, SEED_UNIT_INIT_CAP, inum);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  u->len  = 0;
  u->cap  = SEED_UNIT_INIT_CAP;
  u->inum = inum;
}

void seed_unit_grow(SeedUnit *u, int64 need)
{ int64 newcap = u->cap;
  int64 target = u->len + need;
  if (newcap == 0)
    newcap = need;
  while (newcap < target)
    newcap *= 2;
  uint8 *p = realloc(u->data, newcap);
  if (p == NULL)
    { int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      fprintf(stderr, "[rank %d] seed_unit_grow: realloc(%lld) failed for dest %d\n",
              rank, (long long)newcap, u->inum);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  u->data = p;
  u->cap  = newcap;
}

void seed_unit_trim(SeedUnit *u)
{ if (u->len == 0)
    { free(u->data);
      u->data = NULL;
      u->cap  = 0;
      return;
    }
  uint8 *p = realloc(u->data, u->len);
  if (p != NULL)
    u->data = p;
  u->cap = u->len;
}

void seed_unit_free(SeedUnit *u)
{ free(u->data);
  u->data = NULL;
  u->len  = 0;
  u->cap  = 0;
}
