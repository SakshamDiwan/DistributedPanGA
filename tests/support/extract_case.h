#ifndef PGA_TEST_EXTRACT_CASE_H
#define PGA_TEST_EXTRACT_CASE_H

// Helpers shared by the two extraction binaries.
//
// src/extract.c reads a contig in windows of SCAN_MAX bases, re-reading the
// next chunk when it runs out, with delicate pointer arithmetic to keep
// offsets consistent across the seam.  In production SCAN_MAX is 10 million,
// so no realistic test contig would ever cross a seam.
//
// Hence two binaries: test_extract at the production SCAN_MAX, and
// test_extract_window built with -DSCAN_MAX=64, where a 200-base contig
// crosses several seams.
//
// Both run the SAME fixtures (WINDOW_FIXTURE_SEEDS below) against the SAME
// independent oracle.  That is stronger than diffing the two binaries against
// each other: if both pass, the single-window and multi-window paths
// necessarily agree, AND a bug that corrupted both paths identically still
// could not slip through, because the oracle is independent of both.

#include <stdlib.h>
#include <string.h>
#include "check.h"
#include "synth_gdb.h"
#include "ref_syncmer.h"
#include "rec_build.h"
#include "extract.h"
#include "record.h"

static inline void extract_list(SynthGDB *sg, RefRecList *got)
{
    RecordBuffer buf;
    RecordSizing s = compute_record_sizing(&sg->gdb);
    record_buffer_init(&buf, s);
    extract_to_records(&sg->gdb, 1, 0, &buf);
    rb_decode_all(&buf, got);
    record_buffer_free(&buf);
}

static inline void oracle_list(const char *const *seqs, int nseq, RefRecList *want)
{
    int i;
    ref_list_init(want);
    for (i = 0; i < nseq; i++)
        ref_extract_contig(seqs[i], (int) strlen(seqs[i]), i, want);
}

// Compare two record sets as MULTISETS of (kmer, position, contig, strand):
// sort both, walk in step.  Order-independent on purpose -- extraction emits
// in scan order while the oracle emits in offset order -- and duplicates must
// be preserved, so neither a list compare nor a set compare would do.
static inline void expect_same_multiset(RefRecList *got, RefRecList *want,
                                        const char *what)
{
    int64_t i, lim;
    ref_list_sort(got);
    ref_list_sort(want);
    CHECK_MSG(got->n == want->n, "%s: got %lld records, want %lld",
              what, (long long) got->n, (long long) want->n);
    lim = (got->n < want->n) ? got->n : want->n;
    for (i = 0; i < lim; i++) {
        if (ref_rec_cmp(&got->v[i], &want->v[i]) != 0) {
            t_checks++;
            t_fail_at(__FILE__, __LINE__, "record multiset mismatch",
                      "%s: first diff at %lld: got (c%d pos%lld s%d) want (c%d pos%lld s%d)",
                      what, (long long) i,
                      got->v[i].contig,  (long long) got->v[i].position,  got->v[i].strand,
                      want->v[i].contig, (long long) want->v[i].position, want->v[i].strand);
            return;
        }
    }
    CHECK(1);
}

// Build a fixed-seed contig of the given length, extract it, and require an
// exact match against the oracle.
static inline void extract_case_vs_oracle(int len, uint32_t seed, const char *label)
{
    char       *seq = synth_random_seq(len, seed);
    const char *seqs[1];
    SynthGDB    sg;
    RefRecList  got, want;

    seqs[0] = seq;
    CHECK(synth_gdb_build(&sg, seqs, 1) == 0);
    extract_list(&sg, &got);
    oracle_list(seqs, 1, &want);
    expect_same_multiset(&got, &want, label);
    ref_list_free(&got);
    ref_list_free(&want);
    synth_gdb_free(&sg);
    free(seq);
}

// Fixtures shared by both binaries, so the two scan paths are compared on
// identical input. 200 bases forces >= 3 windows when SCAN_MAX == 64.
#define WINDOW_FIXTURE_LEN 200
static const uint32_t WINDOW_FIXTURE_SEEDS[] = { 11u, 12345u, 99999u, 271828u, 1u };
#define WINDOW_FIXTURE_N \
    ((int) (sizeof WINDOW_FIXTURE_SEEDS / sizeof WINDOW_FIXTURE_SEEDS[0]))

#endif
