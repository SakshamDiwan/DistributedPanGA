// ============================================================================
// P2 -- Local sorting and semantic LCP.  Tests S-01..S-05b, plus R-04.
// See tests/CATALOG.md for the catalog entries.
// ============================================================================
//
// WHAT "LCP" MEANS HERE, because there are three different things with that
// name in this codebase and only one of them is the contract:
//
//   1. SEMANTIC LCP (what we test): byte 0 of each record holds the number of
//      BASES its k-mer shares with the previous record's k-mer, 0..40.
//      Produced by recalc_all_lcps() and used by the downstream FastGA merge.
//   2. msd_sort's internal byte 0: a different encoding that mixes two
//      formulas and disagrees with (1) at bucket boundaries.  recalc_all_lcps
//      deliberately OVERWRITES it.
//   3. GIXshow's "lcp" column: a raw display of (2).  Not a valid oracle.
//
// So: sort_records() alone does NOT produce semantic LCPs -- it leaves (2)
// behind.  Only recalc_all_lcps() does, and in production it is called from
// run_stage2() right after the sort.  Tests here must not confuse the two.
//
// BUILD FLAG: these tests must be compiled with -DLCPs, which the Makefile
// supplies via $(CFLAGS) because LCPS defaults to 1.  Without it
// recalc_all_lcps takes a different branch that writes a 0/1 boundary marker
// instead of a base count, and S-02 fails with 0 and 1 where it expects 40
// and 17.  (That happened during development and looked like a real bug.)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "synth_gdb.h"
#include "ref_syncmer.h"
#include "rec_build.h"

#include "sort.h"
#include "extract.h"
#include "record.h"
#include "constants.h"

static RecordSizing SZ;   // yeast7-shaped: post=3, cont=1, record_size=16

static void mk(uint8_t *k, int byte, int val)
{
    memset(k, 0, KBYTES);
    if (byte >= 0) k[byte] = (uint8_t) val;
}

// S-01 -- simple_kmer_lcp_bases counts matching BASES, not matching bytes.
//
// Each byte holds four bases, two bits each, first base in the high bits.  So
// a difference in one byte can mean an LCP anywhere in a four-base range, and
// the function has to look inside the byte.  The table below walks one case
// per within-byte position plus both ends of the 40-base range:
//
//   bytes equal everywhere                -> 40  (identical k-mers)
//   byte 0 differs in bits 7..6           ->  0  (base 0 differs)
//   byte 0 differs in bits 1..0           ->  3  (bases 0,1,2 match)
//   byte 1 differs in bits 7..6           ->  4  (all of byte 0 matched)
//   byte 1 differs in bits 1..0           ->  7
//   byte 2 differs in bits 7..6           ->  8  (a byte boundary)
//   byte 9 differs in bits 1..0           -> 39  (last possible mismatch)
//
// mk() writes a single byte and zeroes the rest, so e.g. 0x40 = 01 00 00 00
// differs from 0x00 in the FIRST base of that byte, and 0x01 = 00 00 00 01
// differs in the LAST base.
//
// The hand values are then cross-checked twice: against ref_lcp_bases (which
// unpacks both k-mers to base strings and counts a plain character prefix),
// and over 4000 random pairs biased to share leading bytes so long LCPs occur.
static void s01_lcp_table(void)
{
    uint8_t a[KBYTES], b[KBYTES];
    struct { int byte; int aval; int bval; int want; const char *why; } cases[] = {
        { -1, 0,    0,    40, "identical"      },
        {  0, 0x00, 0x40,  0, "differ base 0"  },
        {  0, 0x00, 0x01,  3, "differ base 3"  },
        {  1, 0x00, 0x40,  4, "differ base 4"  },
        {  1, 0x00, 0x01,  7, "differ base 7"  },
        {  2, 0x00, 0x80,  8, "differ base 8"  },
        {  9, 0x00, 0x01, 39, "differ base 39" },
    };
    int i;

    for (i = 0; i < (int) (sizeof cases / sizeof cases[0]); i++) {
        mk(a, cases[i].byte, cases[i].aval);
        mk(b, cases[i].byte, cases[i].bval);
        CHECK_MSG(simple_kmer_lcp_bases(a, b, KBYTES) == cases[i].want,
                  "%s: got %d want %d", cases[i].why,
                  simple_kmer_lcp_bases(a, b, KBYTES), cases[i].want);
        // The independent reference must agree with the hand value too.
        CHECK_MSG(ref_lcp_bases(a, b) == cases[i].want,
                  "reference disagrees for %s", cases[i].why);
    }

    // Randomised cross-check: production vs independent reference.
    {
        uint32_t x = 987654321u;
        int t;
        for (t = 0; t < 4000; t++) {
            int j;
            for (j = 0; j < KBYTES; j++) {
                x ^= x << 13; x ^= x >> 17; x ^= x << 5;
                a[j] = (uint8_t) (x & 0xff);
                x ^= x << 13; x ^= x >> 17; x ^= x << 5;
                // bias toward shared prefixes so long LCPs occur
                b[j] = (j < 4 && (x & 1)) ? a[j] : (uint8_t) (x & 0xff);
            }
            CHECK_MSG(simple_kmer_lcp_bases(a, b, KBYTES) == ref_lcp_bases(a, b),
                      "production/reference LCP disagree at t=%d", t);
        }
    }
}

// S-02 -- recalc_all_lcps fills byte 0 for a whole buffer.
//
// Contract (src/sort.h): overwrite byte 0 of every record with the number of
// bases shared with its PREDECESSOR; the first record has none, so it gets 0.
//
// Fixture: three records, k-mers [k1, k1, k2], where k2 first differs from k1
// at base 17.  Base 17 lives in byte 4 (bases 16..19) at within-byte position
// 1, so k2 is k1 with byte 4 set to 0x10 = 00 01 00 00.
//
// Expected: 0 (no predecessor), 40 (identical k-mers), 17 (first difference).
// The LCP bytes are pre-set to 99 beforehand, so the test also proves the
// function overwrites rather than merely filling in blanks.
static void s02_recalc_multi(void)
{
    RecordBuffer buf;
    uint8_t k1[KBYTES], k2[KBYTES], lcp;
    RefRec r;

    memset(k1, 0, KBYTES);
    memset(k2, 0, KBYTES);
    k2[4] = 0x10;                       // bases 16..19 = 0,1,0,0 -> first differ at 17
    CHECK_EQ_I(ref_lcp_bases(k1, k2), 17);

    rb_alloc_exact(&buf, SZ, 3);
    rb_set(&buf, 0, k1, 100, 0, 0, 99);   // pre-set LCPs must be overwritten
    rb_set(&buf, 1, k1, 200, 0, 0, 99);
    rb_set(&buf, 2, k2, 300, 0, 0, 99);

    recalc_all_lcps(&buf);

    rb_decode(&buf, 0, &r, &lcp); CHECK_EQ_I(lcp, 0);
    rb_decode(&buf, 1, &r, &lcp); CHECK_EQ_I(lcp, 40);
    rb_decode(&buf, 2, &r, &lcp); CHECK_EQ_I(lcp, 17);
    rb_free(&buf);
}

// S-03(a) -- An empty buffer is a no-op.  This half is correct today.
//
// The singleton half of S-03 is a different story and lives in R-04 below:
// count == 0 and count == 1 are currently handled by the same early return,
// but they should not be.
static void s03_recalc_empty(void)
{
    RecordBuffer buf;
    rb_alloc_exact(&buf, SZ, 0);
    recalc_all_lcps(&buf);          // must not crash or write
    CHECK_EQ_I(buf.count, 0);
    rb_free(&buf);
}

// R-04 / S-03(b) -- A one-record buffer must still have its LCP byte set.
//
// THIS TEST IS EXPECTED TO FAIL until the defect is fixed.  It asserts the
// CORRECT behaviour, not the current behaviour.
//
// The contract in src/sort.h is "overwrite byte 0 of EVERY record".  A single
// record has no predecessor, so 0 is the only consistent value -- exactly as
// the first record of a longer buffer gets 0 in S-02.
//
// The defect: src/sort.c:40 reads `if (buf->count <= 1) return;`, which lumps
// together two cases that differ.  For count == 0 returning immediately is
// right.  For count == 1 it leaves byte 0 holding whatever was there before --
// and after sort_records() that is msd_sort's internal encoding, not 0.  The
// fixture pre-sets the byte to 37 to stand in for that leftover value.
//
// Latent today because run_stage2 rarely hands it a single-record slice -- but
// the empty-rank work in R-01 makes exactly that case reachable.
//
// Fix:  if (count == 0) return;  data[0] = 0;  if (count == 1) return;
static void r04_recalc_singleton(void)
{
    RecordBuffer buf;
    uint8_t k[KBYTES], lcp = 0xFF;
    RefRec r;

    memset(k, 0x24, KBYTES);
    rb_alloc_exact(&buf, SZ, 1);
    rb_set(&buf, 0, k, 7, 0, 0, 37);    // stand-in for whatever msd_sort left
    recalc_all_lcps(&buf);
    rb_decode(&buf, 0, &r, &lcp);

    EXPECT_FAIL_UNTIL("R-04", lcp == 0);
    if (lcp != 0)
        printf("      [observed: singleton LCP byte left at %u; "
               "src/sort.c:40 returns early for count <= 1]\n", lcp);
    rb_free(&buf);
}

// S-04 -- sort_records sorts by k-mer and neither loses nor invents records.
//
// Two independent properties, because neither implies the other:
//   ORDERING     adjacent k-mers must be non-decreasing.  A sort that dropped
//                half the records would still satisfy this.
//   CONSERVATION the multiset of (kmer, position, contig, strand) must be
//                unchanged.  A no-op that sorted nothing would satisfy this.
//
// Byte 0 is EXCLUDED from the comparison: msd_sort legitimately rewrites the
// LCP byte as it sorts, so comparing whole records would fail for a reason
// that has nothing to do with conservation.
//
// The fixture draws 5000 records from a pool of only 1000 distinct k-mers, and
// every seventh record reuses a fixed (position, contig, strand), so both
// duplicate KEYS (ties the sort must handle) and duplicate WHOLE RECORDS
// (which a set-based implementation would silently collapse) occur.  Since
// that depends on chance, the test counts them and asserts both are non-zero
// rather than assuming.
static void s04_sort_orders_and_preserves(void)
{
    enum { N = 5000, POOL = 1000 };
    RecordBuffer buf;
    RefRecList before, after;
    uint8_t (*pool)[KBYTES] = malloc(POOL * KBYTES);
    uint32_t x = 24680u;
    int64_t i;

    for (i = 0; i < POOL; i++) {
        int j;
        for (j = 0; j < KBYTES; j++) {
            x ^= x << 13; x ^= x >> 17; x ^= x << 5;
            pool[i][j] = (uint8_t) (x & 0xff);
        }
    }

    rb_alloc_exact(&buf, SZ, N);
    for (i = 0; i < N; i++) {
        int64_t pick;
        x ^= x << 13; x ^= x >> 17; x ^= x << 5;
        pick = x % POOL;
        // i % 7 == 0 repeats the previous occurrence exactly -> duplicate records
        if (i > 0 && (i % 7) == 0)
            rb_set(&buf, i, pool[(x + 1) % POOL], 42, 3, 1, 0);
        else
            rb_set(&buf, i, pool[pick], (int64_t) (i * 13 % 100000),
                   (int) (i % 100), (int) (i & 1), 0);
    }

    rb_decode_all(&buf, &before);

    // Anti-vacuity guard.  The fixture is built from a small pool so that
    // duplicates arise by collision, which is a probabilistic claim -- so
    // check it rather than assume it.  Without duplicate k-mers the ordering
    // check never exercises ties; without duplicate whole records the
    // conservation check never distinguishes a multiset from a set.
    {
        RefRecList probe;
        int64_t j, dup_kmer = 0, dup_rec = 0;
        rb_decode_all(&buf, &probe);
        ref_list_sort(&probe);
        for (j = 1; j < probe.n; j++) {
            if (memcmp(probe.v[j - 1].kmer, probe.v[j].kmer, KBYTES) == 0) dup_kmer++;
            if (ref_rec_cmp(&probe.v[j - 1], &probe.v[j]) == 0)            dup_rec++;
        }
        CHECK_MSG(dup_kmer > 0, "fixture has no duplicate k-mers");
        CHECK_MSG(dup_rec  > 0, "fixture has no duplicate occurrences");
        printf("      [fixture: %lld duplicate-kmer adjacencies, %lld identical occurrences]\n",
               (long long) dup_kmer, (long long) dup_rec);
        ref_list_free(&probe);
    }

    sort_records(&buf);
    rb_decode_all(&buf, &after);

    // ordering: adjacent k-mers non-decreasing
    for (i = 1; i < buf.count; i++) {
        if (memcmp(after.v[i - 1].kmer, after.v[i].kmer, KBYTES) > 0) {
            t_checks++;
            t_fail_at(__FILE__, __LINE__, "sorted order",
                      "k-mer order violated at record %lld", (long long) i);
            break;
        }
    }
    CHECK(1);

    // occurrence multiset preserved (kmer, position, contig, strand)
    ref_list_sort(&before);
    ref_list_sort(&after);
    CHECK_EQ_I(after.n, before.n);
    for (i = 0; i < before.n && i < after.n; i++) {
        if (ref_rec_cmp(&before.v[i], &after.v[i]) != 0) {
            t_checks++;
            t_fail_at(__FILE__, __LINE__, "occurrence multiset",
                      "first diff at %lld", (long long) i);
            break;
        }
    }
    CHECK(1);

    ref_list_free(&before); ref_list_free(&after);
    rb_free(&buf); free(pool);
}

// S-05a -- The sorter writes at most one byte past the last record.
//
// msd_sort writes a sentinel byte at array[nelem * rsize], one past the end of
// the records.  Every buffer that reaches it must therefore be allocated with
// an extra byte.  src/extract.c does that with the `+ 1` in
// record_buffer_init and record_buffer_grow.
//
// This test allocates EXACTLY count*record_size + 1 bytes with no spare
// capacity, so the sentinel write lands on that final byte and nowhere else.
// Under ASan, a sorter that wrote two bytes past, or that wrote past a full
// buffer, would be caught immediately.
//
// SCOPE LIMIT, stated because it matters: the test allocates the +1 ITSELF.
// So it constrains sort_records only.  If someone deleted the +1 from the
// production allocators, this test would still pass.  Catching that is
// S-05b's job, and S-05b was verified to actually catch it.
static void s05a_sorter_sentinel(void)
{
    enum { N = 64 };
    RecordBuffer buf;
    uint32_t x = 5150u;
    int64_t i;

    rb_alloc_exact(&buf, SZ, N);          // exactly N*rsize + 1, no spare capacity
    for (i = 0; i < N; i++) {
        uint8_t k[KBYTES]; int j;
        for (j = 0; j < KBYTES; j++) {
            x ^= x << 13; x ^= x >> 17; x ^= x << 5;
            k[j] = (uint8_t) (x & 0xff);
        }
        rb_set(&buf, i, k, i, 0, 0, 0);
    }
    CHECK_EQ_I(buf.count, buf.capacity);
    sort_records(&buf);                   // ASan catches any write past the +1
    CHECK(1);
    rb_free(&buf);
}

// S-05b -- PRODUCTION allocation and growth really do reserve the sentinel.
//
// Where S-05a tests the sorter, this tests the allocators, and it never
// allocates anything itself: the buffer comes from record_buffer_init()
// followed by extract_to_records(), which is exactly how a real rank builds
// one.
//
// The fixture is 1.5 Mbp so that the record count passes the initial capacity
// of 1<<20 and forces record_buffer_grow() to run at least once -- otherwise
// only the init path would be covered.  (Observed: count grows to ~1.19M and
// capacity doubles to 2^21.)
//
// It then writes to data[capacity * record_size] -- precisely the byte the
// production allocators promise to have reserved.  Under ASan that write is
// a heap-buffer-overflow if either allocator dropped its `+ 1`.
//
// VERIFIED TO FAIL: during development the `+ 1` was deleted from both
// allocators and this test reported a heap-buffer-overflow at exactly this
// line, then passed again once restored.  A test of an invariant like this is
// worth nothing unless you have watched it fail.
//
// Skipped by default (needs --slow) because generating and scanning 1.5 Mbp
// takes a few seconds; `make test-asan` passes --slow.
static void s05b_production_alloc_sentinel(void)
{
    const int LEN = 1500000;              // enough records to force a growth
    char *seq = synth_random_seq(LEN, 777u);
    const char *seqs[1] = { seq };
    SynthGDB sg;
    RecordBuffer buf;
    RecordSizing s;

    CHECK(synth_gdb_build(&sg, seqs, 1) == 0);
    s = compute_record_sizing(&sg.gdb);
    record_buffer_init(&buf, s);
    CHECK_EQ_I(buf.capacity, 1 << 20);    // production's initial capacity
    extract_to_records(&sg.gdb, 1, 0, &buf);

    CHECK_MSG(buf.capacity > (1 << 20),
              "fixture did not force record_buffer_grow (count=%lld cap=%lld)",
              (long long) buf.count, (long long) buf.capacity);
    printf("      [production buffer: count=%lld capacity=%lld rsize=%d]\n",
           (long long) buf.count, (long long) buf.capacity, s.record_size);

    // The sentinel byte production promises: capacity*record_size is in bounds.
    buf.data[buf.capacity * s.record_size] = 0xAA;
    CHECK_EQ_I(buf.data[buf.capacity * s.record_size], 0xAA);

    record_buffer_free(&buf);
    synth_gdb_free(&sg);
    free(seq);
}

int main(int argc, char **argv)
{
    int slow = (argc > 1 && strcmp(argv[1], "--slow") == 0);
    SZ = rs_make(3, 1);
    CHECK_EQ_I(SZ.record_size, 16);

    printf("test_lcp_sort\n");
    RUN(s01_lcp_table);
    RUN(s02_recalc_multi);
    RUN(s03_recalc_empty);
    RUN(r04_recalc_singleton);
    RUN(s04_sort_orders_and_preserves);
    RUN(s05a_sorter_sentinel);
    if (slow) RUN(s05b_production_alloc_sentinel);
    else printf("  %-52s %s\n", "s05b_production_alloc_sentinel", "skipped (--slow)");
    return t_report("test_lcp_sort");
}
