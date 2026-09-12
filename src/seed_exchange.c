#include <limits.h>
#include "seed_exchange.h"

int plan_byte_exchange(const int64_t *sizes, int n, int *counts, int *displs)
{
    int64_t cum = 0;
    int     r;

    for (r = 0; r < n; r++) {
        int64_t s = sizes[r];

        if (s < 0)                    return -1;
        if (s > (int64_t) INT_MAX)    return -2;
        // The running displacement must fit too: every destination can be
        // under INT_MAX while their sum is not, which would wrap displs[]
        // negative and make MPI read outside the send buffer.
        if (cum + s > (int64_t) INT_MAX) return -3;

        counts[r] = (int) s;
        displs[r] = (int) cum;
        cum += s;
    }
    return 0;
}
