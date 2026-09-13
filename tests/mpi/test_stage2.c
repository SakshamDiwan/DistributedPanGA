// ============================================================================
// P5 -- The real distributed sort, under real MPI.  Tests M-01..M-06 and the
// R-01 seam regression.  See tests/CATALOG.md.
// ============================================================================
//
// MUST BE RUN UNDER srun, INSIDE A SLURM ALLOCATION.  tests/run_mpi.sh drives
// it at several world sizes.  Without srun, cray-mpich initialises a singleton
// (world size 1), which still exercises the uniform and skew modes but cannot
// reach any cross-rank behaviour.
//
// WHAT run_stage2 DOES, since everything here is checking one call to it: each
// rank arrives holding an arbitrary slice of records.  After the call, records
// have been redistributed so that rank 0 holds the globally smallest k-mers,
// rank 1 the next, and so on; each rank's slice is locally sorted; and every
// record's LCP byte holds the base-overlap with its global predecessor --
// including across rank boundaries, which requires ranks to exchange their
// edge k-mers.
//
// NO FIXTURES NEEDED.  Records are generated from a PRNG seeded by their
// GLOBAL index, so rank r simply generates the slice it owns.  Two consequences
// make the tests possible:
//   - the expected output multiset is just "records 0..N-1", regenerable
//     anywhere without communication;
//   - the union of all slices is IDENTICAL at every world size, which is what
//     lets M-06 compare a 2-rank run against an 8-rank run at all.
//
// MODES (argv[1]):
//   uniform  N records split evenly.  Covers M-01 (nothing lost or invented),
//            M-02 (global ordering), M-03 (LCP correctness including seams),
//            M-05 (repeatability), and writes the digest M-06 compares.
//   skew     M-04: rank 0 starts with zero records and rank 1 with exactly one.
//   seam     R-01: forces an EMPTY rank between two populated ones.
//            Requires exactly 4 ranks.
//
// EXIT STATUS is Allreduced so every rank returns the same code; otherwise
// srun reports a confusing mix of successes and failures for one logical run.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "check.h"
#include "rec_build.h"
#include "ref_syncmer.h"

#include "dsort.h"
#include "sort.h"
#include "extract.h"
#include "record.h"
#include "constants.h"

static RecordSizing SZ;
static int RANK, WS;

// The k-mer for global record index g.  Deterministic and independent of rank
// or world size, which is what makes the expected output regenerable anywhere
// and makes cross-world-size comparison (M-06) meaningful.
static void gen_kmer(int64_t g, uint8_t *k)
{
    uint32_t x = (uint32_t) (g * 2654435761u) ^ 0x9E3779B9u;
    int j;
    for (j = 0; j < KBYTES; j++) {
        x ^= x << 13; x ^= x >> 17; x ^= x << 5;
        k[j] = (uint8_t) (x & 0xff);
    }
}

static void fill_slice(RecordBuffer *b, int64_t lo, int64_t hi)
{
    int64_t i;
    rb_alloc_exact(b, SZ, hi - lo);
    for (i = lo; i < hi; i++) {
        uint8_t k[KBYTES];
        gen_kmer(i, k);
        rb_set(b, i - lo, k, i % 1000000, (int) (i % 64), (int) (i & 1), 0);
    }
}

// Gather every rank's post-sort slice onto rank 0, concatenated IN RANK ORDER.
//
// Rank order matters: the global-ordering and LCP checks are statements about
// the sequence you get by reading rank 0's slice, then rank 1's, and so on.
// Plain MPI_Gatherv is used rather than production's mpi_gather_records, which
// belongs to the superseded stage-1 path and is not under test.
static uint8_t *gather_all(const RecordBuffer *r, int64_t *n_out)
{
    int      rsize = SZ.record_size;
    int      mine  = (int) (r->count * rsize);
    int     *cnt   = NULL, *dsp = NULL;
    uint8_t *all   = NULL;
    int64_t  total = 0;
    int      i;

    if (RANK == 0) { cnt = malloc(WS * sizeof(int)); dsp = malloc(WS * sizeof(int)); }
    MPI_Gather(&mine, 1, MPI_INT, cnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (RANK == 0) {
        for (i = 0; i < WS; i++) { dsp[i] = (int) total; total += cnt[i]; }
        all = malloc((size_t) (total ? total : 1));
    }
    MPI_Gatherv(r->data, mine, MPI_BYTE, all, cnt, dsp, MPI_BYTE, 0, MPI_COMM_WORLD);
    if (RANK == 0) { *n_out = total / rsize; free(cnt); free(dsp); }
    return all;
}

// Decode a gathered blob into a RefRecList, keeping the LCP bytes alongside.
static void decode_blob(const uint8_t *blob, int64_t n, RefRecList *out, uint8_t **lcps)
{
    RecordBuffer view;
    int64_t i;
    memset(&view, 0, sizeof view);
    view.sizing = SZ; view.count = n; view.capacity = n;
    view.data = (uint8_t *) blob;
    ref_list_init(out);
    out->v = malloc((size_t) (n ? n : 1) * sizeof(RefRec));
    out->cap = n; out->n = n;
    *lcps = malloc((size_t) (n ? n : 1));
    for (i = 0; i < n; i++) rb_decode(&view, i, &out->v[i], &(*lcps)[i]);
}

// Run M-01, M-02 and M-03 against a gathered result.  Rank 0 only.
//
// The three are independent and none implies the others: ordering could hold
// on a truncated result, conservation could hold on an unsorted one, and LCP
// bytes could be wrong on a correctly sorted, complete result.
static void check_global(const uint8_t *blob, int64_t n, int64_t total_in,
                         const char *what)
{
    RefRecList got, want;
    uint8_t   *lcps;
    int64_t    i;

    decode_blob(blob, n, &got, &lcps);

    // M-02: globally non-decreasing k-mer order across the rank-order concatenation
    for (i = 1; i < n; i++)
        if (memcmp(got.v[i - 1].kmer, got.v[i].kmer, KBYTES) > 0) {
            t_checks++;
            t_fail_at(__FILE__, __LINE__, "global k-mer order",
                      "%s: violated at %lld", what, (long long) i);
            break;
        }
    CHECK(1);

    // M-03: every LCP byte must equal the base-overlap with the record before
    // it in the GLOBAL sequence -- which includes the records that sit either
    // side of a rank boundary, the case R-01 is about.
    //
    // Checked with ref_lcp_bases, the independent reference that unpacks both
    // k-mers to base strings and counts a plain character prefix.  Using
    // production's simple_kmer_lcp_bases here would be self-referential: it is
    // the very function run_stage2 used to produce these bytes, so the two
    // would agree even if both were wrong.
    if (n > 0) CHECK_MSG(lcps[0] == 0, "%s: first record LCP = %u, want 0", what, lcps[0]);
    for (i = 1; i < n; i++) {
        int want_lcp = ref_lcp_bases(got.v[i - 1].kmer, got.v[i].kmer);
        if (lcps[i] != want_lcp) {
            t_checks++;
            t_fail_at(__FILE__, __LINE__, "semantic LCP",
                      "%s: record %lld has LCP %u, want %d",
                      what, (long long) i, lcps[i], want_lcp);
            break;
        }
    }
    CHECK(1);

    // M-01: occurrence multiset equals the generated input
    ref_list_init(&want);
    want.v = malloc((size_t) (total_in ? total_in : 1) * sizeof(RefRec));
    want.cap = total_in; want.n = total_in;
    for (i = 0; i < total_in; i++) {
        gen_kmer(i, want.v[i].kmer);
        want.v[i].position = i % 1000000;
        want.v[i].contig   = (int) (i % 64);
        want.v[i].strand   = (int) (i & 1);
    }
    CHECK_MSG(n == total_in, "%s: got %lld records, generated %lld",
              what, (long long) n, (long long) total_in);
    ref_list_sort(&got); ref_list_sort(&want);
    for (i = 0; i < n && i < total_in; i++)
        if (ref_rec_cmp(&got.v[i], &want.v[i]) != 0) {
            t_checks++;
            t_fail_at(__FILE__, __LINE__, "occurrence multiset",
                      "%s: first diff at %lld", what, (long long) i);
            break;
        }
    CHECK(1);

    ref_list_free(&got); ref_list_free(&want); free(lcps);
}

// Write a canonical digest of the result so run_mpi.sh can compare across
// world sizes (M-06).
//
// "Canonical" means the SORTED occurrence multiset.  It deliberately does not
// record the order in which equal k-mers appear: records sharing a k-mer may
// be permuted among themselves depending on which rank contributed them, and
// stable ordering of equal keys is not a documented contract of the sort.
// Requiring byte-identical output would make this test fail for a legal
// implementation detail.
static void write_digest(const uint8_t *blob, int64_t n, const char *path)
{
    RefRecList l; uint8_t *lcps; FILE *f; int64_t i;
    decode_blob(blob, n, &l, &lcps);
    ref_list_sort(&l);                    // canonical: sorted occurrence multiset
    f = fopen(path, "w");
    if (f == NULL) { perror("digest"); ref_list_free(&l); free(lcps); return; }
    fprintf(f, "records %lld\n", (long long) n);
    for (i = 0; i < n; i++) {
        int j;
        for (j = 0; j < KBYTES; j++) fprintf(f, "%02x", l.v[i].kmer[j]);
        fprintf(f, " %lld %d %d\n", (long long) l.v[i].position,
                l.v[i].contig, l.v[i].strand);
    }
    fclose(f);
    ref_list_free(&l); free(lcps);
}

// ---- modes -----------------------------------------------------------------

// Mode "uniform" -- the baseline: N records spread evenly over the ranks.
static void mode_uniform(const char *digest)
{
    const int64_t TOTAL = 4000;
    RecordBuffer local, recv;
    uint8_t *blob; int64_t n = 0;
    int64_t lo = (TOTAL * RANK) / WS, hi = (TOTAL * (RANK + 1)) / WS;

    fill_slice(&local, lo, hi);
    memset(&recv, 0, sizeof recv);
    run_stage2(&local, &recv, MPI_COMM_WORLD);
    rb_free(&local);

    CHECK_EQ_I(verify_per_rank_sorted(&recv, NULL), 0);   // M-02, per-rank half

    blob = gather_all(&recv, &n);
    if (RANK == 0) {
        check_global(blob, n, TOTAL, "M-01..M-03 uniform");
        if (digest) write_digest(blob, n, digest);
        free(blob);
    }
    record_buffer_free(&recv);

    // M-05: repeatability.  Run the identical input a second time and require
    // the same SEMANTIC invariants -- conservation, ordering, LCP values.
    //
    // Deliberately not required: byte-identical slices.  See write_digest
    // above; equal k-mers have no guaranteed relative order, so demanding
    // byte equality would pin an accident of the sort rather than a contract.
    {
        RecordBuffer l2, r2;
        uint8_t *b2; int64_t n2 = 0;
        fill_slice(&l2, lo, hi);
        memset(&r2, 0, sizeof r2);
        run_stage2(&l2, &r2, MPI_COMM_WORLD);
        rb_free(&l2);
        b2 = gather_all(&r2, &n2);
        if (RANK == 0) { check_global(b2, n2, TOTAL, "M-05 repeat"); free(b2); }
        record_buffer_free(&r2);
    }
}

// Mode "skew" -- M-04: degenerate slice sizes on the INPUT side.
//
// Rank 0 starts with zero records and rank 1 with exactly one; the rest get a
// normal share.  Empty and single-element buffers are the classic place for
// off-by-one and null-pointer faults in redistribution code.
//
// Note the asymmetry being tested: a rank that starts empty may still RECEIVE
// records, and a rank that starts full may end up empty, because ownership is
// decided by k-mer value and not by who supplied the data.  So all the
// assertions are about the post-exchange result, not the input.
static void mode_skew(void)
{
    const int64_t PER = 700;
    RecordBuffer local, recv;
    uint8_t *blob; int64_t n = 0, lo, hi, total = 0;
    int r;

    // rank 0 -> 0 records, rank 1 -> 1 record, others -> PER
    for (r = 0; r < WS; r++) total += (r == 0) ? 0 : (r == 1 ? 1 : PER);
    lo = 0;
    for (r = 0; r < RANK; r++) lo += (r == 0) ? 0 : (r == 1 ? 1 : PER);
    hi = lo + ((RANK == 0) ? 0 : (RANK == 1 ? 1 : PER));

    fill_slice(&local, lo, hi);
    CHECK_EQ_I(local.count, hi - lo);
    memset(&recv, 0, sizeof recv);
    run_stage2(&local, &recv, MPI_COMM_WORLD);
    rb_free(&local);

    CHECK_EQ_I(verify_per_rank_sorted(&recv, NULL), 0);
    blob = gather_all(&recv, &n);
    if (RANK == 0) { check_global(blob, n, total, "M-04 skew"); free(blob); }
    record_buffer_free(&recv);
}

// Mode "seam" -- R-01: an EMPTY rank sitting between two populated ranks.
//
// FIXED.  This test was landed as an expected failure, reproduced the defect
// under real MPI, and the marker was removed only after the fix made it XPASS.
//
// THE DEFECT WAS.  After sorting, each rank needs the last k-mer of the rank
// before it, to compute the LCP of its own first record.  run_stage2 obtained
// it with a single MPI_Sendrecv between NEIGHBOURS -- it only ever looked one
// hop.  A rank holding no records has no last k-mer, so it sent an all-zero
// sentinel.  Its successor then computed an LCP against 0x00... instead of
// against the real predecessor, sitting two or more ranks back.  The boundary
// CHECK passed vacuously, so nothing aborted and a wrong LCP byte flowed
// downstream into seed matching.
//
// Observed before the fix: rank 2 reported LCP 0 where 4 is correct.  The same
// defect also surfaced through mode_skew at world size 2, where the single
// global record lands on rank 1 with rank 0 empty.
//
// Fixed by replacing the one-hop exchange with an MPI_Allgather of every rank's
// (has-records, last-k-mer), from which each rank selects its nearest NON-EMPTY
// predecessor (src/dsort.c).
//
// BUILDING A FIXTURE THAT FORCES AN EMPTY MIDDLE RANK.  Ownership is decided
// by bucket, and a bucket is the first 5 bases: (kmer[0] << 2) | (kmer[1] >> 6).
// Give every record kmer[0] = 0xC0 and put mass in only two buckets:
//
//     kmer[1] >> 6 == 0  ->  bucket 0xC0*4 + 0 = 768
//     kmer[1] >> 6 == 2  ->  bucket 0xC0*4 + 2 = 770
//
// leaving bucket 769 with nothing in it.  With 4 ranks and equal counts,
// compute_destination_map assigns 768 -> rank 0, 769 -> rank 1, 770 -> rank 2.
// Rank 1 therefore owns only the empty bucket: an empty rank between two
// populated ones.  (test_dsort.c S-08 asserts that same ownership as a pure
// function, so the fixture is verified independently of MPI.)
//
// WHY THE WRONG ANSWER AND THE RIGHT ANSWER DIFFER HERE.  Both groups share
// kmer[0] = 0xC0, i.e. their first four bases are identical, and they diverge
// in the top bits of kmer[1] -- base 4.  So the correct LCP is 4.  The buggy
// value compares against all-zero bytes, which differ from 0xC0 at base 0,
// giving 0.  Had the two groups differed in kmer[0], both answers would have
// been 0 and the bug would have hidden.
//
// PRECONDITIONS ARE ASSERTED FIRST -- the ownership map, the slice sizes, and
// that rank 1 is genuinely empty between two populated ranks.  If the splitter
// ever changes, this test must fail loudly as "precondition not met" rather
// than quietly passing because the setup stopped reproducing the bug.
static void mode_seam(void)
{
    enum { PER = 4 };
    RecordBuffer local, recv;
    int64_t counts[8];
    int      i;

    if (WS != 4) {
        if (RANK == 0) printf("  %-52s %s\n", "r01_seam", "skipped (needs -n 4)");
        return;
    }

    // Ranks 0 and 1 supply one bucket group each, and ranks 2 and 3 supply
    // nothing.  Who ends up empty is then decided purely by the ownership map,
    // not by which rank happened to produce the data -- so the test is really
    // about redistribution rather than about the input layout.
    rb_alloc_exact(&local, SZ, (RANK < 2) ? PER : 0);
    for (i = 0; i < ((RANK < 2) ? PER : 0); i++) {
        uint8_t k[KBYTES];
        memset(k, 0, KBYTES);
        k[0] = 0xC0;
        k[1] = (uint8_t) ((RANK == 0 ? 0x00 : 0x80) | i);   // bucket 768 or 770
        k[2] = (uint8_t) (0x10 + i);                         // distinct within bucket
        k[9] = (uint8_t) i;
        rb_set(&local, i, k, 100 + i, RANK, 0, 0);
    }

    memset(&recv, 0, sizeof recv);
    run_stage2(&local, &recv, MPI_COMM_WORLD);
    rb_free(&local);

    // Precondition 1: ownership is what the fixture assumes.
    if (RANK == 0) {
        int64_t hist[NUM_BUCK];
        int    *ks = malloc(5 * sizeof(int)), *sel = malloc(NUM_BUCK * sizeof(int));
        memset(hist, 0, sizeof hist);
        hist[768] = PER; hist[770] = PER;
        compute_destination_map(hist, 4, ks, sel);
        CHECK_MSG(sel[768] == 0, "precondition: select[768] = %d, want 0", sel[768]);
        CHECK_MSG(sel[769] == 1, "precondition: select[769] = %d, want 1", sel[769]);
        CHECK_MSG(sel[770] == 2, "precondition: select[770] = %d, want 2", sel[770]);
        free(ks); free(sel);
    }

    // Precondition 2: rank 1 really is empty, between two populated ranks.
    {
        int64_t mine = recv.count;
        MPI_Gather(&mine, 1, MPI_LONG_LONG, counts, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        if (RANK == 0) {
            printf("      [slice sizes: r0=%lld r1=%lld r2=%lld r3=%lld]\n",
                   (long long) counts[0], (long long) counts[1],
                   (long long) counts[2], (long long) counts[3]);
            CHECK_MSG(counts[0] == PER, "precondition: rank 0 should hold %d", PER);
            CHECK_MSG(counts[1] == 0,   "precondition: rank 1 should be empty");
            CHECK_MSG(counts[2] == PER, "precondition: rank 2 should hold %d", PER);
        }
    }

    // The assertion itself.  Rank 0 sends its LAST k-mer directly to rank 2,
    // bypassing the empty rank, and rank 2 compares its own first record's LCP
    // byte against the correct answer computed from those two k-mers.
    //
    // Correct value: 4 (shared kmer[0] = 0xC0, diverging at base 4).
    // Buggy value:   0 (compared against rank 1's all-zero sentinel).
    {
        uint8_t last0[KBYTES], first2[KBYTES], lcp2 = 0;
        MPI_Status st;

        if (RANK == 0 && recv.count > 0)
            MPI_Send(recv.data + (recv.count - 1) * SZ.record_size + KMER_OFFSET,
                     KBYTES, MPI_BYTE, 2, 77, MPI_COMM_WORLD);
        if (RANK == 2) {
            MPI_Recv(last0, KBYTES, MPI_BYTE, 0, 77, MPI_COMM_WORLD, &st);
            memcpy(first2, recv.data + KMER_OFFSET, KBYTES);
            lcp2 = recv.data[LCP_OFFSET];
            {
                int want = ref_lcp_bases(last0, first2);
                printf("      [rank2 seam: observed LCP %u, correct %d]\n", lcp2, want);
                CHECK_MSG(lcp2 == want,
                          "seam LCP past an empty rank is %u, want %d", lcp2, want);
            }
        }
    }
    record_buffer_free(&recv);
}

int main(int argc, char **argv)
{
    const char *mode = (argc > 1) ? argv[1] : "uniform";
    const char *digest = (argc > 2) ? argv[2] : NULL;
    int local_fail, global_fail;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &WS);
    SZ = rs_make(3, 1);

    if (RANK == 0) printf("test_stage2 mode=%s world=%d\n", mode, WS);

    if      (strcmp(mode, "uniform") == 0) mode_uniform(digest);
    else if (strcmp(mode, "skew")    == 0) mode_skew();
    else if (strcmp(mode, "seam")    == 0) mode_seam();
    else { if (RANK == 0) printf("unknown mode %s\n", mode); MPI_Finalize(); return 2; }

    // Every rank must exit with the same status, or srun reports a confusing mix.
    local_fail = t_fails;
    MPI_Allreduce(&local_fail, &global_fail, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (RANK == 0) {
        printf("test_stage2[%s,ws=%d]: %d checks, %d failures", mode, WS,
               t_checks, global_fail - t_xpass);
        if (t_xfail) printf(", %d xfail", t_xfail);
        printf("\n");
    }
    MPI_Finalize();
    return global_fail ? 1 : 0;
}
