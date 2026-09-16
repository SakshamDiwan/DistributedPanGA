// ============================================================================
// P1 -- Record field widths and whole-contig work splitting.
// Tests E-12, E-13, and the R-03 sizing table.  See tests/CATALOG.md.
// ============================================================================
//
// TWO UNRELATED THINGS LIVE HERE because both are pure arithmetic over GDB
// metadata and neither needs sequence data:
//
//   compute_record_sizing (src/record.c)  -- how many bytes each variable-width
//       field of a record gets, derived from the largest contig and the contig
//       count.  Get this wrong and every record silently truncates a value.
//
//   compute_work_split (src/work_split.c) -- which contigs each rank extracts
//       from.  Proposal item O2b generalises this to split WITHIN a contig, so
//       these invariants are what that change must preserve.
//
// Links only record.c and work_split.c: both read GDB struct fields but call
// no GDB functions, so the GDB structs can be filled in by hand and no lib/
// sources are needed.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "record.h"
#include "work_split.h"
#include "constants.h"

static RecordSizing sizing_for(int64_t maxctg, int ncontig)
{
    GDB g;
    memset(&g, 0, sizeof g);
    g.maxctg  = maxctg;
    g.ncontig = ncontig;
    return compute_record_sizing(&g);
}

// How many bytes are needed to represent every value in [0, max_value]?
//
// Written from first principles -- b bytes hold 0 .. 256^b - 1 -- specifically
// so it does NOT share logic with production's bytes_needed().  The whole
// point of R-03 is that the two disagree, so deriving the expectation from the
// code under test would pin the bug instead of catching it.
static int bytes_for_max_value(int64_t max_value)
{
    int     b   = 0;
    int64_t cap = 0;                       // 256^0 - 1 == 0
    while (cap < max_value) { b++; cap = (cap + 1) * 256 - 1; }
    return b ? b : 1;                      // a 0-byte field can represent nothing
}

// R-03 (width half) -- post_bytes must be able to hold every position that
//                       extraction can actually emit.
//
// PART OF THIS TEST IS EXPECTED TO FAIL until the defect is fixed.  The table
// states the CORRECT widths, not the current ones.
//
// The reasoning.  It is tempting to assume positions run 0..maxctg-1, which is
// what bytes_needed(maxctg) computes.  But a reverse-complement record stores
// j + TMER, and j itself can reach len - TMER, so:
//
//     largest stored position = (len - 12) + 12 = len = maxctg
//
// Positions therefore span [0, maxctg] INCLUSIVE -- maxctg + 1 distinct
// values -- and the field must satisfy maxctg <= 256^post_bytes - 1.
//
// The off-by-one only bites when maxctg is an exact power of 256, because only
// then does the extra value cross a byte boundary.  At maxctg = 256 the field
// gets 1 byte, which holds 0..255, and position 256 wraps to 0.
//
// test_extract.c contains the other half of R-03: a real 256-base contig whose
// emitted position 256 is observed decoding as 0.
static void r03_post_bytes_representability(void)
{
    struct { int64_t maxctg; int intended; } cases[] = {
        {     1, 1 }, {   255, 1 }, {   256, 2 }, {   257, 2 },
        { 65535, 2 }, { 65536, 3 }, { 65537, 3 },
    };
    int i;

    for (i = 0; i < (int) (sizeof cases / sizeof cases[0]); i++) {
        int64_t m   = cases[i].maxctg;
        int     got = sizing_for(m, 8).post_bytes;
        int     want = cases[i].intended;

        CHECK_MSG(bytes_for_max_value(m) == want,
                  "table self-check failed for maxctg=%lld", (long long) m);

        if (got == want) {
            CHECK_MSG(got == want, "maxctg=%lld", (long long) m);
        } else {
            // Defect bites exactly at maxctg = 256^k, where bytes_needed's
            // `while (cum < max_val)` stops one byte short.
            printf("      [maxctg=%lld: post_bytes=%d, need %d to hold position %lld]\n",
                   (long long) m, got, want, (long long) m);
            EXPECT_FAIL_UNTIL("R-03", got == want);
        }
    }
}

// The neighbouring field is CORRECT, and this pins that so the R-03 fix does
// not "helpfully" change it too.
//
// cont_bytes = bytes_needed(2 * ncontig).  The factor of 2 looks arbitrary but
// is exactly right: the strand flag occupies the field's top bit, so the ids
// only have 8*b - 1 bits available, and doubling the count reserves that bit.
//   ncontig = 128 -> 2*128 = 256 -> 1 byte: ids 0..127 in bits 0..6, strand
//                    in bit 7.  The tightest legal fit.
//   ncontig = 129 -> 2*129 = 258 -> 2 bytes, correctly.
static void r03_cont_bytes_ok(void)
{
    // ids 0..ncontig-1 must fit below the strand bit at 8*b - 1.
    CHECK_EQ_I(sizing_for(1000, 127).cont_bytes, 1);
    CHECK_EQ_I(sizing_for(1000, 128).cont_bytes, 1);   // ids 0..127 + bit 7
    CHECK_EQ_I(sizing_for(1000, 129).cont_bytes, 2);
    CHECK_EQ_I(sizing_for(1000,   1).cont_bytes, 1);
}

// E-13 -- record_size is the sum of its parts, for the shape real data takes.
//
// Layout (src/record.h): 1 LCP byte + 10 k-mer bytes + 1 mask byte, then the
// two variable-width fields.  For a yeast-sized genome (largest contig ~1.5 M
// bases, ~120 contigs) that gives post_bytes 3 and cont_bytes 1, hence
// record_size 16 -- the shape most of the other tests use as their fixture.
static void e13_record_size(void)
{
    RecordSizing s = sizing_for(1500000, 120);
    CHECK_EQ_I(s.post_bytes, 3);
    CHECK_EQ_I(s.cont_bytes, 1);
    CHECK_EQ_I(s.record_size, 1 + KBYTES + 1 + 3 + 1);
    CHECK_EQ_I(s.record_size, 16);
}

// ---- E-12: whole-contig work split ----------------------------------------

// E-12 helper -- the invariants a work split must satisfy for any input.
//
// compute_work_split divides contigs among workers by cumulative base count,
// at whole-contig boundaries only.  Deliberately NOT asserted: how evenly it
// balances.  That is a policy O2b will change.  What must survive is that the
// split is a genuine PARTITION -- every contig extracted by exactly one rank.
// Extract a contig twice and seeds are duplicated; miss one and alignments
// silently disappear.
static void check_split(const int64_t *clen, int ncontig, int nw, const char *what)
{
    GDB      g;
    int     *split = malloc((size_t) (nw + 1) * sizeof(int));
    int64_t *post  = malloc((size_t) (nw + 1) * sizeof(int64_t));
    int64_t  total = 0;
    int      i, w;

    memset(&g, 0, sizeof g);
    g.ncontig = ncontig;
    g.contigs = calloc((size_t) ncontig + 1, sizeof(GDB_CONTIG));
    for (i = 0; i < ncontig; i++) { g.contigs[i].clen = clen[i]; total += clen[i]; }
    g.seqtot = total;

    compute_work_split(&g, nw, split, post);

    CHECK_MSG(split[0] == 0, "%s: split[0] != 0", what);
    CHECK_MSG(split[nw] == ncontig, "%s: split[nw]=%d != ncontig=%d",
              what, split[nw], ncontig);
    CHECK_MSG(post[0] == 0, "%s: post[0] != 0", what);
    CHECK_MSG(post[nw] == total, "%s: post[nw]=%lld != seqtot=%lld",
              what, (long long) post[nw], (long long) total);

    for (w = 0; w < nw; w++)
        CHECK_MSG(split[w] <= split[w + 1],
                  "%s: split not non-decreasing at %d (%d > %d)",
                  what, w, split[w], split[w + 1]);

    // Every contig owned exactly once.
    {
        int *owned = calloc((size_t) ncontig, sizeof(int));
        for (w = 0; w < nw; w++)
            for (i = split[w]; i < split[w + 1]; i++)
                if (i >= 0 && i < ncontig) owned[i]++;
        for (i = 0; i < ncontig; i++)
            CHECK_MSG(owned[i] == 1, "%s: contig %d owned %d times", what, i, owned[i]);
        free(owned);
    }

    free(g.contigs); free(split); free(post);
}

// E-12 -- Whole-contig work splitting across awkward inputs.
//
// The four shapes are chosen to break different assumptions:
//   (a) eight equal contigs    -- the easy case, and with 3 workers a count
//                                 that does not divide evenly
//   (b) one 10 Mbp contig plus 20 tiny ones -- the balance is hopeless (a
//                                 single contig cannot be shared), so this
//                                 checks the partition still holds when the
//                                 policy cannot do anything sensible.  This is
//                                 precisely the case O2b exists to fix.
//   (c) a single contig        -- including with 4 workers, where three of
//                                 them must legitimately get nothing
//   (d) more workers than contigs -- empty ranges must be representable rather
//                                 than an error or an out-of-range index
static void e12_work_split(void)
{
    {   // (a) eight equal contigs
        int64_t c[8]; int i;
        for (i = 0; i < 8; i++) c[i] = 1000;
        check_split(c, 8, 4, "8 equal / 4 workers");
        check_split(c, 8, 3, "8 equal / 3 workers");
    }
    {   // (b) one dominant contig plus many small ones
        int64_t c[21]; int i;
        c[0] = 10000000;
        for (i = 1; i < 21; i++) c[i] = 500;
        check_split(c, 21, 4, "1 giant + 20 small / 4 workers");
    }
    {   // (c) a single contig
        int64_t c[1] = { 5000 };
        check_split(c, 1, 1, "1 contig / 1 worker");
        check_split(c, 1, 4, "1 contig / 4 workers");   // empty ranges are legitimate
    }
    {   // (d) more workers than contigs
        int64_t c[3] = { 100, 200, 300 };
        check_split(c, 3, 8, "3 contigs / 8 workers");
    }
}

int main(void)
{
    printf("test_record_split\n");
    RUN(r03_post_bytes_representability);
    RUN(r03_cont_bytes_ok);
    RUN(e13_record_size);
    RUN(e12_work_split);
    return t_report("test_record_split");
}
