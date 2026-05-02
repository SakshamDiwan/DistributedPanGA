#include "work_split.h"

void compute_work_split(GDB *gdb, int num_workers, int *split, int64_t *post)
{
    int64_t total_bases    = gdb->seqtot;
    int64_t bases_so_far   = 0;
    int64_t next_threshold = total_bases / num_workers;
    int     next_worker    = 1;
    int     r;

    split[0] = 0;
    post[0]  = 0;

    for (r = 0; r < gdb->ncontig; r++)
    {
        bases_so_far += gdb->contigs[r].clen;
        while (bases_so_far >= next_threshold && next_worker < num_workers)
        {
            split[next_worker] = r + 1;
            post[next_worker]  = bases_so_far;
            next_worker += 1;
            next_threshold = (total_bases * next_worker) / num_workers;
        }
    }
    split[num_workers] = gdb->ncontig;
    post[num_workers]  = total_bases;
}
