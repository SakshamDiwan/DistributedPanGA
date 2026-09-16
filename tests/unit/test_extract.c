// ============================================================================
// P1 -- Extraction, coordinates, strand.  Tests E-01..E-11, plus the R-03
// extraction regression.  See tests/CATALOG.md for the catalog entries.
// ============================================================================
//
// HOW THESE TESTS ARE PUT TOGETHER
//
//   1. A synthetic GDB is built in memory from an ACGT string
//      (tests/support/synth_gdb.c).  No files on disk, no FastGA tools.
//   2. Production extract_to_records() runs on it, exactly as a real rank
//      would.
//   3. The resulting binary records are decoded by an INDEPENDENT decoder
//      (tests/support/rec_build.c), written from the layout spec in
//      src/record.h rather than by calling any production decode path.
//   4. The decoded records are compared against an INDEPENDENT oracle
//      (tests/support/ref_syncmer.c) that re-derives which k-mers should have
//      been emitted, from the algorithm definition.
//
// Why go through extract_to_records() rather than the scanner directly:
// scan_contig(), record_emit_cb() and the EmitFn typedef are all `static` /
// file-local in src/extract.c, so a test cannot call or intercept them.  The
// public entry point is extract_to_records(), so that is what we drive.
//
// Expected values are always derived by hand or by the oracle -- never
// captured from a run of the code under test.  A test that records current
// behaviour would happily pin a bug in place.  That is exactly what R-03
// below demonstrates: the code and the contract disagree, and the test sides
// with the contract.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "synth_gdb.h"
#include "ref_syncmer.h"
#include "rec_build.h"
#include "extract_case.h"

#include "extract.h"
#include "record.h"
#include "pack.h"
#include "tables.h"
#include "constants.h"

// ---------------------------------------------------------------- helpers --

// Search for a pseudorandom sequence in which the 12-mer at offset j_target is
// a selected syncmer.  Used by fixtures that must have a record at one exact
// offset (E-09, R-03).  Roughly 2 in 5 offsets are selected, so this succeeds
// almost immediately; the bounded loop just avoids an infinite hang.  The
// caller asserts the resulting precondition rather than trusting it.
static char *seq_with_selected_j(int len, int64_t j_target, uint32_t *seed_out)
{
    uint32_t seed;
    for (seed = 1; seed < 20000u; seed++) {
        char *s = synth_random_seq(len, seed);
        if (s != NULL && ref_is_selected(s, len, j_target)) {
            if (seed_out) *seed_out = seed;
            return s;
        }
        free(s);
    }
    return NULL;
}

// ------------------------------------------------------------------ tests --

// E-01 -- Forward 40-mer packing is 2 bits per base, most significant first.
//
// Contract (src/record.h, src/pack.c): a 40-mer occupies 10 bytes, four bases
// per byte, with the FIRST base in the HIGH bits of each byte:
//
//     byte = (base0 << 6) | (base1 << 4) | (base2 << 2) | base3
//
// Fixture: "ACGT" repeated 13 times (52 bases).  Every aligned group of four
// bases is therefore A,C,G,T.  With the encoding A=0 C=1 G=2 T=3 that packs to
//
//     00 01 10 11  (binary)  =  0x1B
//
// in all ten bytes.  52 bases is just a convenient length; the 40-mer at
// offset 0 only needs 40.
static void e01_forward_packing(void)
{
    char seq[53];
    const char *seqs[1] = { seq };
    SynthGDB sg;
    uint8_t out[KBYTES];
    int i;

    for (i = 0; i < 13; i++) memcpy(seq + 4 * i, "ACGT", 4);
    seq[52] = '\0';

    CHECK(synth_gdb_build(&sg, seqs, 1) == 0);
    CHECK(pack_forward_kmer(&sg.gdb, 0, 0, out) == 0);
    // A=0,C=1,G=2,T=3 -> 00 01 10 11b = 0x1B
    for (i = 0; i < KBYTES; i++) CHECK_EQ_I(out[i], 0x1B);
    synth_gdb_free(&sg);
}

// E-02 -- Reverse-complement packing both complements and reverses.
//
// Contract (src/pack.c): pack_rc_kmer(gdb, c, j, out) reads the FORWARD bases
// in [j-28, j+12) and stores their reverse complement.  The window is offset
// like that because the RC k-mer is keyed to the syncmer at j, not to j-28.
//
// Fixture: a contig of 40 'A's.  Its reverse complement is 40 'T's, and T=3,
// so every packed byte is 11 11 11 11 = 0xFF.
//
// Two independent ways to see the same answer, both used here:
//   - biologically: RC(AAAA...) = TTTT... -> 0xFF per byte;
//   - from the table: init_comp_table() (src/tables.c) builds Comp[] so that
//     Comp[0x00] = 192|48|12|3 = 0xFF, i.e. the RC of four A's is four T's.
//
// Why 40 bases exactly: RC needs j >= KMER-TMER = 28 and j+12 <= clen, so
// j = 28 requires clen >= 40.  40 is the shortest contig that can emit an RC
// record at all.
static void e02_rc_packing_all_a(void)
{
    char seq[41];
    const char *seqs[1] = { seq };
    SynthGDB sg;
    uint8_t out[KBYTES];
    int i;

    memset(seq, 'A', 40);
    seq[40] = '\0';

    CHECK(synth_gdb_build(&sg, seqs, 1) == 0);
    CHECK(pack_rc_kmer(&sg.gdb, 0, 28, out) == 0);
    // RC(AAAA...) = TTTT...; independently Comp[0x00] = 192|48|12|3 = 0xFF.
    for (i = 0; i < KBYTES; i++) CHECK_EQ_I(out[i], 0xFF);
    synth_gdb_free(&sg);
}

// E-03 -- Packing matches an independent packer on non-symmetric input.
//
// E-01 and E-02 use inputs whose answers are uniform (every byte identical),
// which would not catch a byte-order or within-byte-order mistake.  This test
// uses pseudorandom sequence, where any such mistake changes the result, and
// compares against ref_pack_forward / ref_pack_rc in tests/support.
//
// The reference packers are written from the layout definition: for RC it
// builds the complement string explicitly (RC[i] = complement of
// seq[j + 11 - i]) and then packs it, rather than reusing production's
// Comp[] lookup.  So agreement here is genuine cross-checking, not two calls
// into the same code.
static void e03_rc_packing_asymmetric(void)
{
    char *seq = synth_random_seq(100, 12345u);
    const char *seqs[1] = { seq };
    SynthGDB sg;
    uint8_t got[KBYTES], want[KBYTES];
    int64_t j;

    CHECK(synth_gdb_build(&sg, seqs, 1) == 0);
    for (j = 28; j <= 88; j += 7) {
        CHECK(pack_rc_kmer(&sg.gdb, 0, j, got) == 0);
        ref_pack_rc(seq, j, want);
        CHECK_MSG(memcmp(got, want, KBYTES) == 0, "RC mismatch at j=%lld", (long long) j);
    }
    // and forward, for symmetry of coverage
    for (j = 0; j <= 60; j += 7) {
        CHECK(pack_forward_kmer(&sg.gdb, 0, j, got) == 0);
        ref_pack_forward(seq, j, want);
        CHECK_MSG(memcmp(got, want, KBYTES) == 0, "fwd mismatch at j=%lld", (long long) j);
    }
    synth_gdb_free(&sg);
    free(seq);
}

// E-04 -- Forward packing refuses a window that runs off the contig end.
//
// pack_forward_kmer needs the whole 40-mer [j, j+40) inside the contig.
// On a 52-base contig:
//     j = 12 -> needs [12, 52)  -- fits exactly, must succeed
//     j = 13 -> needs [13, 53)  -- one base past the end, must return -1
// Pinning both sides of the boundary is the point; testing only the failing
// side would pass even if the function rejected everything.
static void e04_forward_range(void)
{
    char *seq = synth_random_seq(52, 7u);
    const char *seqs[1] = { seq };
    SynthGDB sg;
    uint8_t out[KBYTES];

    CHECK(synth_gdb_build(&sg, seqs, 1) == 0);
    CHECK(pack_forward_kmer(&sg.gdb, 0, 12, out) == 0);   // [12,52) fits exactly
    CHECK(pack_forward_kmer(&sg.gdb, 0, 13, out) == -1);  // [13,53) overruns
    synth_gdb_free(&sg);
    free(seq);
}

// E-05 -- RC packing refuses windows off EITHER end of the contig.
//
// pack_rc_kmer reads forward bases [j-28, j+12), so it has two ways to fall
// off a 52-base contig:
//     j = 27 -> begb = -1   (runs off the start)   must return -1
//     j = 28 -> begb =  0   first valid offset     must succeed
//     j = 40 -> endb = 52   last valid offset      must succeed
//     j = 41 -> endb = 53   (runs off the end)     must return -1
// Both boundaries are pinned from both sides.
static void e05_rc_range(void)
{
    char *seq = synth_random_seq(52, 9u);
    const char *seqs[1] = { seq };
    SynthGDB sg;
    uint8_t out[KBYTES];

    CHECK(synth_gdb_build(&sg, seqs, 1) == 0);
    CHECK(pack_rc_kmer(&sg.gdb, 0, 27, out) == -1);  // begb = -1
    CHECK(pack_rc_kmer(&sg.gdb, 0, 28, out) ==  0);
    CHECK(pack_rc_kmer(&sg.gdb, 0, 40, out) ==  0);
    CHECK(pack_rc_kmer(&sg.gdb, 0, 41, out) == -1);  // endb = 53
    synth_gdb_free(&sg);
    free(seq);
}

// E-06 -- Extraction matches an independent syncmer oracle.  THE O1 GATE.
//
// This is the most important test in the suite.  Proposal item O1 replaces the
// extraction path with HySortK; this test is what will detect any change to
// which k-mers get emitted, at which positions, on which strand.
//
// The contract, restated from the algorithm rather than from src/extract.c:
// FastGA selects CLOSED (12,8)-SYNCMERS.  For the 12-mer starting at offset j,
// look at its five constituent 8-mers, at offsets j, j+1, j+2, j+3, j+4.  Hash
// each one canonically (the smaller of the forward and reverse-complement
// hash, so the choice does not depend on which strand you read).  The 12-mer
// is selected if and only if the SMALLEST of those five hashes sits at either
// END of the window -- position j or position j+4.  Ties select.
//
// Each selected j can emit up to two records:
//     forward, if j <= len-40 : stored position j,    strand 0
//     RC,      if j >= 28     : stored position j+12, strand 1
//
// tests/support/ref_syncmer.c implements that definition naively and
// separately.  It shares only TMap[] -- a fixed scramble table that is part of
// the algorithm's definition, like a magic constant -- and recomputes the
// 4-base reverse complement itself instead of reusing production's Comp[].
//
// ANTI-VACUITY: the production scanner reaches its selection decision through
// three different branches (a new strict minimum at the right end; the old
// minimum falling off the left end, forcing a rescan; and an exact tie).  If a
// fixture never exercised one of them, this test could pass while that branch
// was broken.  So the oracle counts how many offsets hit each branch and we
// assert all three are non-zero.  On this fixture the tie branch fires only
// twice in 4000 bases -- which is precisely why counting is worth doing.
static void e06_oracle_match(void)
{
    char *seq = synth_random_seq(4000, 20260910u);
    const char *seqs[1] = { seq };
    SynthGDB sg;
    RefRecList got, want;

    CHECK(synth_gdb_build(&sg, seqs, 1) == 0);
    extract_list(&sg, &got);
    oracle_list(seqs, 1, &want);

    // Anti-vacuity: every selection branch must actually be reached, or a
    // pass would prove nothing about the tie / rescan paths.
    CHECK_MSG(want.n_right > 0, "right-end branch never exercised");
    CHECK_MSG(want.n_left  > 0, "left-end rescan branch never exercised");
    CHECK_MSG(want.n_tie   > 0, "tie fall-through never exercised");
    printf("      [oracle: %lld selected, right=%lld left=%lld tie=%lld, %lld records]\n",
           (long long) want.n_sel, (long long) want.n_right,
           (long long) want.n_left, (long long) want.n_tie, (long long) want.n);

    expect_same_multiset(&got, &want, "E-06 4kb contig");

    ref_list_free(&got); ref_list_free(&want);
    synth_gdb_free(&sg); free(seq);
}

// E-07 -- The two strands use different position conventions.
//
// This is separated from E-06 so that a coordinate bug reports as a
// coordinate failure rather than as an opaque "record sets differ".
//
// Contract (src/extract.c, the two emit() calls):
//     forward record: stored position is j itself
//     RC record:      stored position is j + TMER  (= j + 12)
//
// The asymmetry is deliberate -- the RC k-mer covers [j-28, j+12), so its
// natural coordinate is the right end of that window, not the left.
//
// The check runs in reverse: for each decoded record, recover the offset j it
// claims to come from (j = position for forward, position - 12 for RC) and ask
// the oracle whether that offset really is a selected syncmer.  So it also
// catches records invented at offsets that should not have been selected.
static void e07_coordinates(void)
{
    char *seq = synth_random_seq(2000, 4242u);
    const char *seqs[1] = { seq };
    SynthGDB sg;
    RefRecList got;
    int64_t i, nf = 0, nr = 0;
    int len = 2000;

    CHECK(synth_gdb_build(&sg, seqs, 1) == 0);
    extract_list(&sg, &got);

    for (i = 0; i < got.n; i++) {
        RefRec *r = &got.v[i];
        if (r->strand == 0) {
            // forward: position is the 12-mer offset j itself, and must be selected
            CHECK_MSG(ref_is_selected(seq, len, r->position),
                      "fwd record at pos %lld is not a selected offset",
                      (long long) r->position);
            CHECK_MSG(r->position <= len - KMER, "fwd pos %lld > len-40",
                      (long long) r->position);
            nf++;
        } else {
            // RC: position is j + TMER, so j = position - 12 must be selected
            int64_t j = r->position - TMER;
            CHECK_MSG(ref_is_selected(seq, len, j),
                      "rc record at pos %lld implies unselected j=%lld",
                      (long long) r->position, (long long) j);
            CHECK_MSG(j >= KMER - TMER, "rc j %lld < 28", (long long) j);
            nr++;
        }
    }
    CHECK(nf > 0); CHECK(nr > 0);
    ref_list_free(&got); synth_gdb_free(&sg); free(seq);
}

// E-08 -- Near the ends of a short contig, only one strand can be emitted.
//
// The two emission conditions have opposite senses, so on a short contig they
// carve out disjoint ranges of j.  On a 41-base contig:
//     forward needs  j <= len-40 = 1     -> only j = 0 or 1
//     RC needs       j >= 28             -> only j = 28 or 29
//                                           (j <= len-12 = 29)
// Since 1 < 28 the two ranges cannot overlap, so NO offset emits both strands.
// Contrast E-09, which finds the exact contig length at which they can.
static void e08_end_asymmetry(void)
{
    char *seq = synth_random_seq(41, 31337u);
    const char *seqs[1] = { seq };
    SynthGDB sg;
    RefRecList got, want;
    int64_t i;

    CHECK(synth_gdb_build(&sg, seqs, 1) == 0);
    extract_list(&sg, &got);
    oracle_list(seqs, 1, &want);

    for (i = 0; i < got.n; i++) {
        if (got.v[i].strand == 0)
            CHECK_MSG(got.v[i].position <= 1, "fwd j=%lld exceeds len-40=1",
                      (long long) got.v[i].position);
        else
            CHECK_MSG(got.v[i].position - TMER >= 28, "rc j=%lld below 28",
                      (long long) (got.v[i].position - TMER));
    }
    expect_same_multiset(&got, &want, "E-08 41-base contig");
    ref_list_free(&got); ref_list_free(&want);
    synth_gdb_free(&sg); free(seq);
}

// E-09 -- Both strands at the SAME offset become possible at exactly len 68.
//
// For one offset j to emit both records it must satisfy both conditions at
// once:  28 <= j  and  j <= len - 40.  Such a j exists only when
// len - 40 >= 28, i.e. len >= 68.
//
// So this test pins the boundary from both sides, using j = 28:
//     len = 67 -> len-40 = 27, and 28 > 27, so forward is impossible:
//                 RC record at position 40 only
//     len = 68 -> len-40 = 28, so j = 28 just qualifies:
//                 BOTH a forward record at 28 and an RC record at 40
//
// The fixtures are found by searching seeds until the oracle reports that
// offset 28 is a selected syncmer -- otherwise there would be no record at
// that offset at all and the test would pass vacuously.
static void e09_both_orientations_boundary(void)
{
    uint32_t seed67 = 0, seed68 = 0;
    char *s67 = seq_with_selected_j(67, 28, &seed67);
    char *s68 = seq_with_selected_j(68, 28, &seed68);
    int k;

    CHECK_MSG(s67 != NULL && s68 != NULL, "could not build fixtures with j=28 selected");
    if (s67 == NULL || s68 == NULL) { free(s67); free(s68); return; }

    for (k = 0; k < 2; k++) {
        const char *seq = k ? s68 : s67;
        int         len = k ? 68  : 67;
        const char *seqs[1] = { seq };
        SynthGDB sg;
        RefRecList got;
        int64_t i, fwd28 = 0, rc40 = 0;

        CHECK(synth_gdb_build(&sg, seqs, 1) == 0);
        extract_list(&sg, &got);
        for (i = 0; i < got.n; i++) {
            if (got.v[i].strand == 0 && got.v[i].position == 28) fwd28++;
            if (got.v[i].strand == 1 && got.v[i].position == 40) rc40++;
        }
        CHECK_MSG(rc40 == 1, "len=%d: expected one RC record at pos 40, got %lld",
                  len, (long long) rc40);
        if (len == 67)
            CHECK_MSG(fwd28 == 0, "len=67: forward at j=28 must be impossible (28 > 27)");
        else
            CHECK_MSG(fwd28 == 1, "len=68: forward at j=28 must be present");
        ref_list_free(&got); synth_gdb_free(&sg);
    }
    free(s67); free(s68);
}

// E-10 -- Strand rides in the top bit of the contig field without corrupting
//          the contig id, even at the largest id that field can hold.
//
// Contract (src/record.h): the contig field is `cont_bytes` bytes, stored
// little-endian, and the strand flag is packed into the HIGH bit of that
// field:  value = contig_id | (strand << (8*cont_bytes - 1)).
//
// Fixture: 128 contigs.  compute_record_sizing derives cont_bytes from
// 2*ncontig = 256, giving cont_bytes = 1.  That is the tightest legal fit:
// ids 0..127 occupy bits 0..6 and the strand flag takes bit 7.  Contig 127 is
// therefore the maximum id, and the one where a sign-extension or masking
// mistake would show up -- decoding it must give exactly (127, strand), not
// 255 or -1.
static void e10_strand_roundtrip(void)
{
    enum { NC = 128, LEN = 80 };
    char       *seq[NC];
    const char *seqs[NC];
    SynthGDB    sg;
    RecordSizing s;
    RefRecList  got;
    int64_t     i, from127 = 0, s0 = 0, s1 = 0;
    int         c;

    for (c = 0; c < NC; c++) { seq[c] = synth_random_seq(LEN, 1000u + c); seqs[c] = seq[c]; }

    CHECK(synth_gdb_build(&sg, seqs, NC) == 0);
    s = compute_record_sizing(&sg.gdb);
    // 2*ncontig = 256 -> 1 byte: ids 0..127 in bits 0..6, strand in bit 7.
    CHECK_EQ_I(s.cont_bytes, 1);

    extract_list(&sg, &got);
    for (i = 0; i < got.n; i++) {
        CHECK_MSG(got.v[i].contig >= 0 && got.v[i].contig < NC,
                  "decoded contig %d out of range", got.v[i].contig);
        if (got.v[i].contig == 127) {
            from127++;
            if (got.v[i].strand) s1++; else s0++;
        }
    }
    CHECK_MSG(from127 > 0, "no records decoded from contig 127");
    CHECK_MSG(s0 > 0 && s1 > 0,
              "contig 127 (max id) needs both strands: fwd=%lld rc=%lld",
              (long long) s0, (long long) s1);

    ref_list_free(&got); synth_gdb_free(&sg);
    for (c = 0; c < NC; c++) free(seq[c]);
}

// R-03 (extraction half) -- A real extraction that loses a position.
//
// THIS TEST IS EXPECTED TO FAIL until the defect is fixed.  It asserts the
// CORRECT behaviour, not the current behaviour.
//
// The defect: positions are stored in `post_bytes` bytes, and
// compute_record_sizing sets post_bytes = bytes_needed(maxctg).  But the
// largest position actually emitted is maxctg, not maxctg-1, because an RC
// record stores j + TMER and j can reach len - TMER:
//
//     max stored position = (len - 12) + 12 = len = maxctg
//
// So positions span [0, maxctg] -- that is maxctg + 1 distinct values -- and
// one byte only holds 0..255.  At maxctg = 256 the field is one byte short and
// position 256 wraps to 0.  (See test_record_split.c for the same defect
// expressed as a width table.)
//
// Choosing the fixture: we need a contig where some record really does land on
// position 256.  Position 256 requires j + 12 = 256, so j = 244; and j is
// capped at len - 12 = 244.  So j = 244 is the ONLY offset that can produce
// it, and the contig must be exactly 256 bases.  We search seeds until the
// oracle says offset 244 is a selected syncmer, then assert that precondition
// explicitly -- otherwise no record would be emitted there and the test would
// pass while proving nothing.
static void r03_position_256_extraction(void)
{
    uint32_t seed = 0;
    char *seq = seq_with_selected_j(256, 244, &seed);
    const char *seqs[1] = { seq };
    SynthGDB sg;
    RecordSizing s;
    RefRecList got;
    int64_t i, at256 = 0, at0 = 0;

    CHECK_MSG(seq != NULL, "could not build a 256-base fixture with j=244 selected");
    if (seq == NULL) return;

    CHECK(synth_gdb_build(&sg, seqs, 1) == 0);
    s = compute_record_sizing(&sg.gdb);
    CHECK_EQ_I(sg.gdb.maxctg, 256);
    printf("      [fixture seed=%u, maxctg=256, post_bytes=%d]\n", seed, s.post_bytes);

    extract_list(&sg, &got);
    for (i = 0; i < got.n; i++) {
        if (got.v[i].strand == 1 && got.v[i].position == 256) at256++;
        if (got.v[i].strand == 1 && got.v[i].position == 0)   at0++;
    }
    // Precondition: the oracle says this record must exist.
    CHECK_MSG(ref_is_selected(seq, 256, 244), "precondition: j=244 must be selected");

    EXPECT_FAIL_UNTIL("R-03", at256 == 1);
    if (at256 != 1)
        printf("      [observed: RC records decoded at position 0: %lld "
               "(position 256 truncated by post_bytes=%d)]\n",
               (long long) at0, s.post_bytes);

    ref_list_free(&got); synth_gdb_free(&sg); free(seq);
}

// E-11 (default-window half) -- the control group for the window test.
//
// src/extract.c scans a contig in windows of SCAN_MAX bases (10 million in
// production).  Anything longer is re-read in chunks, using some delicate
// pointer arithmetic to keep offsets consistent across the seam.
//
// test_extract_window is built with -DSCAN_MAX=64 so that a 200-base contig
// crosses several seams.  This function runs the SAME fixtures here at the
// production SCAN_MAX, where they fit in a single window.  Both binaries
// compare against the same independent oracle, so if both pass, the
// single-window and multi-window paths necessarily agree with each other --
// and, unlike a direct binary-to-binary diff, a bug that corrupted both paths
// identically still could not slip through.
static void e11_window_fixtures_default(void)
{
    int i;
    for (i = 0; i < WINDOW_FIXTURE_N; i++) {
        char label[64];
        snprintf(label, sizeof label, "E-11 default SCAN_MAX seed=%u",
                 WINDOW_FIXTURE_SEEDS[i]);
        extract_case_vs_oracle(WINDOW_FIXTURE_LEN, WINDOW_FIXTURE_SEEDS[i], label);
    }
}

int main(void)
{
    init_comp_table();
    printf("test_extract\n");
    RUN(e01_forward_packing);
    RUN(e02_rc_packing_all_a);
    RUN(e03_rc_packing_asymmetric);
    RUN(e04_forward_range);
    RUN(e05_rc_range);
    RUN(e06_oracle_match);
    RUN(e07_coordinates);
    RUN(e08_end_asymmetry);
    RUN(e09_both_orientations_boundary);
    RUN(e10_strand_roundtrip);
    RUN(e11_window_fixtures_default);
    RUN(r03_position_256_extraction);
    return t_report("test_extract");
}
