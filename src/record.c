#include "record.h"

// Smallest byte width that can represent every value in [0, max_value].
// Returns b such that 256^b >= max_value + 1, i.e. the field can hold
// max_value itself and not merely max_value - 1.
static int bytes_needed_inclusive(int64_t max_value)
{
    int     bytes = 0;
    int64_t cum   = 1;                    // 256^0 == 1 value: just 0
    while (cum < max_value + 1)
    {
        cum *= 256;
        bytes += 1;
    }
    return bytes;
}

RecordSizing compute_record_sizing(GDB *gdb)
{
    RecordSizing s;

    // Positions run [0, maxctg] INCLUSIVE, not [0, maxctg). A
    // reverse-complement record stores j + TMER and j reaches len - TMER
    // (src/extract.c), so the largest stored position is exactly the contig
    // length. Sizing for maxctg alone is one value short, which only shows up
    // when maxctg is an exact power of 256 -- there a 256-base contig emits
    // position 256 into a one-byte field and it wraps to 0.
    s.post_bytes  = bytes_needed_inclusive(gdb->maxctg);

    // Contig ids run [0, ncontig-1], but the strand flag occupies the field's
    // top bit, leaving 8*b - 1 bits for the id. Doubling the count reserves
    // that bit: 256^b >= 2*ncontig  <=>  2^(8b-1) >= ncontig. Already correct,
    // so it is expressed as an inclusive bound on the largest encoded value.
    s.cont_bytes  = bytes_needed_inclusive((int64_t) 2 * gdb->ncontig - 1);
    s.record_size = 1 + KBYTES + 1 + s.post_bytes + s.cont_bytes;
    return s;
}
