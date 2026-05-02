#include "record.h"

static int bytes_needed(int64_t max_val)
{
    int     bytes = 0;
    int64_t cum   = 1;
    while (cum < max_val)
    {
        cum *= 256;
        bytes += 1;
    }
    return bytes;
}

RecordSizing compute_record_sizing(GDB *gdb)
{
    RecordSizing s;
    s.post_bytes  = bytes_needed(gdb->maxctg);
    s.cont_bytes  = bytes_needed((int64_t) 2 * gdb->ncontig);
    s.record_size = 1 + KBYTES + 1 + s.post_bytes + s.cont_bytes;
    return s;
}
