// ============================================================================
// P2 (distributed half) -- the splitter and the send-buffer packing.
// Tests S-06..S-09.  See tests/CATALOG.md.
// ============================================================================
//
// BACKGROUND.  run_stage2() distributes records across ranks so that the
// global result is sorted.  It does this WITHOUT any rank knowing the whole
// dataset, in four steps:
//
//   1. Each record is assigned a BUCKET from its first 5 bases -- 1024 of
//      them (NUM_BUCK).  Records with the same leading bases share a bucket.
//   2. Every rank counts its own records per bucket, and an MPI_Allreduce
//      turns those into one global histogram that every rank holds.
//   3. Every rank independently runs compute_destination_map() on that same
//      histogram to decide which rank owns which buckets.  Because the input
//      is identical everywhere and the function is deterministic, all ranks
//      reach the same answer with NO further communication.
//   4. organize_send_buffer() groups the local records by destination so a
//      single MPI_Alltoallv can ship them.
//
// Two invariants make the whole scheme work, and both are tested below:
//   - bucket ownership must be MONOTONE (rank 0 owns the lowest buckets, and
//     so on), or the concatenated result would not be globally sorted;
//   - the map must be DETERMINISTIC, or ranks would disagree about who owns
//     what and records would be lost.
//
// These functions are pure.  dsort.h includes mpi.h for its prototypes, so the
// binary is built with mpicc, but nothing here calls MPI: the only MPI call in
// organize_send_buffer is an MPI_Abort on an overflow path we never trigger.
// So this binary runs as an ordinary single process, no srun needed.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "rec_build.h"
#include "ref_syncmer.h"

#include "dsort.h"
#include "extract.h"
#include "record.h"
#include "constants.h"

static RecordSizing SZ;

// S-06 -- The bucket of a record is exactly its first 5 bases, nothing else.
//
// record_bucket (src/dsort.h) computes  (kmer[0] << 2) | (kmer[1] >> 6).
// kmer[0] holds bases 0..3 and the top two bits of kmer[1] hold base 4, so
// that expression is precisely "the first five bases", giving
// 4^5 = 1024 = NUM_BUCK possible values.
//
// Three things are checked:
//   - the formula itself, over all 256 values of kmer[0] and all 4 top-bit
//     patterns of kmer[1];
//   - that the low 6 bits of kmer[1] (bases 5..7) never change the answer --
//     varied explicitly, since using them by accident would still look
//     plausible but would break the monotonicity the splitter relies on;
//   - that the 1024 combinations map onto 1024 DISTINCT buckets covering the
//     whole range, i.e. it is a bijection with no collisions or gaps.
static void s06_record_bucket(void)
{
    uint8_t k[KBYTES];
    int seen[NUM_BUCK];
    int b0, hi, low;

    memset(seen, 0, sizeof seen);
    memset(k, 0, KBYTES);

    for (b0 = 0; b0 < 256; b0++) {
        for (hi = 0; hi < 4; hi++) {
            // low 6 bits of kmer[1] must not matter; vary them
            for (low = 0; low < 64; low += 21) {
                int got;
                k[0] = (uint8_t) b0;
                k[1] = (uint8_t) ((hi << 6) | low);
                got  = record_bucket(k);
                CHECK_MSG(got == ((b0 << 2) | hi),
                          "bucket(%02x,%02x) = %d, want %d", b0, k[1], got, (b0 << 2) | hi);
            }
            k[1] = (uint8_t) (hi << 6);
            seen[record_bucket(k)]++;
        }
    }
    for (b0 = 0; b0 < NUM_BUCK; b0++)
        CHECK_MSG(seen[b0] == 1, "bucket %d hit %d times (want 1)", b0, seen[b0]);
}

// S-07 -- bucket_local_counts tallies a buffer into the per-bucket histogram.
//
// This histogram is what step 2 above feeds into the Allreduce, so an
// off-by-one here would corrupt every rank's view of the data.
//
// The fixture plants known counts in three buckets chosen to cover the range:
// bucket 0 (kmer[0]=0x00, top bits 0), bucket 1023 (kmer[0]=0xFF, top bits 3
// -- the maximum, where an overflow would show), and bucket 512.  Bucket 7 is
// deliberately left empty to check that untouched buckets read back as 0
// rather than as uninitialised memory.  The total is also checked against
// buf->count so nothing is double-counted.
static void s07_bucket_counts(void)
{
    RecordBuffer buf;
    int64_t counts[NUM_BUCK];
    int64_t i, sum = 0;
    // bucket 0 x3, bucket 1023 x2, bucket 512 x1, bucket 7 empty
    struct { int b0, hi; int n; } spec[] = { {0,0,3}, {255,3,2}, {128,0,1} };
    int s, r = 0, total = 0;

    for (s = 0; s < 3; s++) total += spec[s].n;
    rb_alloc_exact(&buf, SZ, total);
    for (s = 0; s < 3; s++) {
        int j;
        for (j = 0; j < spec[s].n; j++) {
            uint8_t k[KBYTES];
            memset(k, 0, KBYTES);
            k[0] = (uint8_t) spec[s].b0;
            k[1] = (uint8_t) (spec[s].hi << 6);
            k[9] = (uint8_t) j;                 // distinct records
            rb_set(&buf, r++, k, j, 0, 0, 0);
        }
    }

    bucket_local_counts(&buf, counts);
    CHECK_EQ_I(counts[0], 3);
    CHECK_EQ_I(counts[(255 << 2) | 3], 2);      // 1023
    CHECK_EQ_I(counts[(128 << 2) | 0], 1);      // bucket 512
    CHECK_EQ_I(counts[7], 0);                   // untouched bucket
    for (i = 0; i < NUM_BUCK; i++) sum += counts[i];
    CHECK_EQ_I(sum, buf.count);
    rb_free(&buf);
}

// S-08 helper -- assert the invariants that must hold for ANY histogram and
// world size.  Called from s08_destination_map with several shapes of input.
//
// Deliberately NOT asserted here: which specific buckets land on which rank.
// That is a balancing policy and O2c is expected to change it.  What must
// survive any policy change are the structural properties below.
static void check_map(const int64_t *hist, int ws, const char *what)
{
    int *ksplit = malloc((size_t) (ws + 1) * sizeof(int));
    int *select = malloc(NUM_BUCK * sizeof(int));
    int *k2     = malloc((size_t) (ws + 1) * sizeof(int));
    int *s2     = malloc(NUM_BUCK * sizeof(int));
    int  b, r;

    compute_destination_map(hist, ws, ksplit, select);
    compute_destination_map(hist, ws, k2, s2);

    // determinism: this is what licenses every rank computing it redundantly
    CHECK_MSG(memcmp(select, s2, NUM_BUCK * sizeof(int)) == 0,
              "%s: select not deterministic", what);
    CHECK_MSG(memcmp(ksplit, k2, (size_t) (ws + 1) * sizeof(int)) == 0,
              "%s: ksplit not deterministic", what);

    for (b = 0; b < NUM_BUCK; b++)
        CHECK_MSG(select[b] >= 0 && select[b] < ws,
                  "%s: select[%d] = %d out of range [0,%d)", what, b, select[b], ws);

    for (b = 1; b < NUM_BUCK; b++)
        CHECK_MSG(select[b - 1] <= select[b],
                  "%s: select not monotone at %d (%d > %d)",
                  what, b, select[b - 1], select[b]);

    CHECK_MSG(ksplit[0] == 0, "%s: ksplit[0] = %d", what, ksplit[0]);
    CHECK_MSG(ksplit[ws] == NUM_BUCK, "%s: ksplit[%d] = %d, want %d",
              what, ws, ksplit[ws], NUM_BUCK);
    for (r = 1; r <= ws; r++)
        CHECK_MSG(ksplit[r - 1] <= ksplit[r],
                  "%s: ksplit not non-decreasing at %d (%d > %d)",
                  what, r, ksplit[r - 1], ksplit[r]);

    free(ksplit); free(select); free(k2); free(s2);
}

// S-08 -- compute_destination_map's structural invariants.  THE O2c GATE.
//
// Proposal item O2c deepens or replaces this splitter to balance better at
// high rank counts.  These are the properties it must not break:
//
//   MONOTONE      select[] non-decreasing.  This is the load-bearing one:
//                 rank r must own a contiguous, ordered range of buckets, or
//                 concatenating the ranks' outputs in rank order would not
//                 produce a globally sorted sequence.
//   IN RANGE      every destination is a real rank in [0, world_size).
//   DETERMINISTIC identical input gives identical output.  This is what
//                 licenses every rank computing the map redundantly instead
//                 of broadcasting it, so losing it would be a silent
//                 correctness bug, not a performance one.
//   CONSISTENT    ksplit[] brackets the same partition, and
//                 ksplit[world_size] == NUM_BUCK.
//
// Input shapes cover the awkward cases: a uniform histogram; all mass in one
// bucket; an entirely EMPTY histogram (legal -- ranks may own nothing, but the
// destination must still be a valid rank); and the two-bucket shape R-01 uses.
// World sizes include 3, because rounding behaviour that looks fine at powers
// of two often does not survive an odd divisor, and 1024 = NUM_BUCK, where
// there is exactly one bucket per rank.
//
// One concrete expectation is pinned: with a perfectly uniform histogram and 4
// ranks, each rank should own exactly a quarter of the buckets, so
// select[b] == b / 256.
static void s08_destination_map(void)
{
    int64_t *hist = calloc(NUM_BUCK, sizeof(int64_t));
    int      b;

    // (a) uniform
    for (b = 0; b < NUM_BUCK; b++) hist[b] = 1;
    check_map(hist, 1, "uniform/1");
    check_map(hist, 3, "uniform/3");
    check_map(hist, 4, "uniform/4");
    check_map(hist, 1024, "uniform/1024");
    {   // each rank owns a contiguous quarter
        int *ks = malloc(5 * sizeof(int)), *sel = malloc(NUM_BUCK * sizeof(int));
        compute_destination_map(hist, 4, ks, sel);
        for (b = 0; b < NUM_BUCK; b++)
            CHECK_MSG(sel[b] == b / 256, "uniform/4: select[%d] = %d, want %d",
                      b, sel[b], b / 256);
        free(ks); free(sel);
    }

    // (b) all mass in one bucket
    memset(hist, 0, NUM_BUCK * sizeof(int64_t));
    hist[500] = 1000;
    check_map(hist, 4, "single bucket/4");
    check_map(hist, 3, "single bucket/3");

    // (c) empty histogram -- empty ranks are legitimate, out-of-range is not
    memset(hist, 0, NUM_BUCK * sizeof(int64_t));
    check_map(hist, 4, "empty/4");
    check_map(hist, 1, "empty/1");

    // (d) two populated buckets, as R-01's fixture uses
    memset(hist, 0, NUM_BUCK * sizeof(int64_t));
    hist[768] = 4; hist[770] = 4;
    check_map(hist, 4, "two buckets/4");
    {   // the exact ownership R-01 depends on
        int *ks = malloc(5 * sizeof(int)), *sel = malloc(NUM_BUCK * sizeof(int));
        compute_destination_map(hist, 4, ks, sel);
        printf("      [R-01 fixture: select[768]=%d select[769]=%d select[770]=%d]\n",
               sel[768], sel[769], sel[770]);
        CHECK_EQ_I(sel[768], 0);
        CHECK_EQ_I(sel[769], 1);
        CHECK_EQ_I(sel[770], 2);
        free(ks); free(sel);
    }
    free(hist);
}

// S-09 -- organize_send_buffer regroups records by destination losslessly.
//
// MPI_Alltoallv needs all records destined for rank r to sit contiguously, at
// a known offset, with a known length.  organize_send_buffer rewrites the
// local buffer into that shape and reports counts[] and displs[] in BYTES.
//
// Checked here:
//   - counts[r] matches an independently computed per-destination tally;
//   - displs[] are the exact prefix sums of counts[], so the regions abut with
//     no gaps or overlaps, and the last one ends at the total size;
//   - every record sits in the slice belonging to ITS destination;
//   - every input record appears exactly once in the output, byte-identical --
//     checked by matching each output record to an unclaimed input record, so
//     both losses and duplications are caught.
//
// Destination 2 is deliberately given no records: an empty destination must
// still get count 0 and a well-defined displacement, since MPI will read those
// entries regardless.
static void s09_organize_send(void)
{
    enum { N = 1000, WS = 3 };
    RecordBuffer buf;
    int     *select = malloc(NUM_BUCK * sizeof(int));
    int      counts[WS], displs[WS];
    uint8_t *send = NULL;
    uint32_t x = 13579u;
    int64_t  i;
    int      b, rsize;
    int64_t  per[WS] = { 0, 0, 0 };

    // destination 2 deliberately receives nothing
    for (b = 0; b < NUM_BUCK; b++) select[b] = (b < 500) ? 0 : 1;

    rb_alloc_exact(&buf, SZ, N);
    rsize = buf.sizing.record_size;
    for (i = 0; i < N; i++) {
        uint8_t k[KBYTES];
        int j;
        for (j = 0; j < KBYTES; j++) {
            x ^= x << 13; x ^= x >> 17; x ^= x << 5;
            k[j] = (uint8_t) (x & 0xff);
        }
        rb_set(&buf, i, k, i, (int) (i % 50), (int) (i & 1), 0);
        per[select[record_bucket(k)]]++;
    }

    organize_send_buffer(&buf, select, WS, &send, counts, displs);

    CHECK_EQ_I(displs[0], 0);
    for (b = 0; b < WS; b++)
        CHECK_MSG(counts[b] == (int) (per[b] * rsize),
                  "dest %d: count %d, want %lld", b, counts[b], (long long) (per[b] * rsize));
    for (b = 1; b < WS; b++)
        CHECK_MSG(displs[b] == displs[b - 1] + counts[b - 1],
                  "dest %d: displ %d != %d + %d", b, displs[b], displs[b - 1], counts[b - 1]);
    CHECK_MSG(counts[2] == 0, "destination 2 should be empty, got %d", counts[2]);
    CHECK_EQ_I(displs[WS - 1] + counts[WS - 1], (int) (N * rsize));

    // every record lands in its destination's slice, exactly once, byte-identical
    {
        int *found = calloc(N, sizeof(int));
        for (b = 0; b < WS; b++) {
            int64_t off;
            for (off = 0; off < counts[b]; off += rsize) {
                const uint8_t *rec = send + displs[b] + off;
                int64_t j;
                int hit = -1;
                CHECK_MSG(select[record_bucket(rec + KMER_OFFSET)] == b,
                          "record in slice %d belongs to %d", b,
                          select[record_bucket(rec + KMER_OFFSET)]);
                for (j = 0; j < N; j++)
                    if (!found[j] && memcmp(rec, buf.data + j * rsize, rsize) == 0) { hit = 1; found[j] = 1; break; }
                if (hit < 0) {
                    t_checks++;
                    t_fail_at(__FILE__, __LINE__, "record not found in input",
                              "slice %d offset %lld", b, (long long) off);
                    break;
                }
            }
        }
        for (i = 0; i < N; i++)
            if (!found[i]) { t_checks++; t_fail_at(__FILE__, __LINE__, "record dropped",
                                                   "input record %lld", (long long) i); break; }
        CHECK(1);
        free(found);
    }

    free(send); free(select); rb_free(&buf);
}

int main(void)
{
    SZ = rs_make(3, 1);
    printf("test_dsort\n");
    RUN(s06_record_bucket);
    RUN(s07_bucket_counts);
    RUN(s08_destination_map);
    RUN(s09_organize_send);
    return t_report("test_dsort");
}
