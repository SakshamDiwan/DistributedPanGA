// ============================================================================
// P4 -- The FastK/FastGA adapter.  Tests K-01..K-05.  See tests/CATALOG.md.
// ============================================================================
//
// WHAT THE ADAPTER IS FOR.  The downstream merge and alignment code is lifted
// from FastGA and expects to read k-mers through FastK's Kmer_Stream interface.
// Our pipeline instead produces a flat array of our own records.  Rather than
// rewrite FastGA, src/kmer_adapter.c builds a KmerStreamAdapter that serves
// the same stream operations (the pipeline uses it directly, no casts).
//
// The compiler checks the struct's fields, but not their MEANING: get the
// index semantics, a width, or the entry layout subtly wrong, and FastGA reads
// plausible garbage instead of failing loudly.  That is what these tests check.
//
// HOW THE STREAM IS ADDRESSED.  A 40-mer is 10 bytes.  The adapter splits it:
//   - the first prefix_bytes = 3 bytes are the PREFIX.  These are not stored
//     per entry at all; the position in cumulative_count_by_prefix implies them.
//   - the remaining entry_kmer_suffix_bytes = 7 are stored with each entry,
//     along with a mask byte, the LCP byte, and the position/contig payload.
// cumulative_count_by_prefix[p] (CBP[p] below) is CUMULATIVE: the number of
// records whose prefix is <= p.  So the entries for prefix p occupy
// [CBP[p-1], CBP[p]), and a reader tracks the current prefix in
// `current_prefix` while walking entries.  Current_Entry() rebuilds the full
// k-mer by prepending current_prefix's three bytes to the stored suffix.
//
// SIZE WARNING: num_prefixes = 2^24 entries of int64 = 128 MiB per adapter, allocated
// regardless of how few records there are.  So each test builds exactly one
// and frees it before the next.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "rec_build.h"
#include "ref_syncmer.h"

#include "kmer_adapter.h"
#include "extract.h"
#include "record.h"
#include "constants.h"

static RecordSizing SZ;

static int prefix_of(const uint8_t *kmer)
{
    return (kmer[0] << 16) | (kmer[1] << 8) | kmer[2];
}

static int cmp_rec(const void *a, const void *b)
{
    return memcmp(((const RefRec *) a)->kmer, ((const RefRec *) b)->kmer, KBYTES);
}

// Build a k-mer-sorted RecordBuffer, and return the same records as a sorted
// array for the test to assert against.
//
// The adapter requires its input already sorted by k-mer (that is what
// run_stage2 hands it), so we sort here with our own comparator rather than
// calling production's sorter -- keeping the model independent.
//
// Every fifth record copies the previous record's 3-byte prefix, which forces
// RUNS of equal prefixes.  Without them the per-prefix counts would be nearly
// all 0s and 1s, largest_prefix_group_num_entries would be trivially 1, and
// K-01 would not distinguish a cumulative count
// from a per-prefix count.  The remaining random prefixes leave large GAPS,
// which exercise the "no records with this prefix" case.
static RefRec *build_sorted(RecordBuffer *buf, int64_t n, uint32_t seed)
{
    RefRec  *m = malloc((size_t) n * sizeof(RefRec));
    uint32_t x = seed;
    int64_t  i;

    for (i = 0; i < n; i++) {
        int j;
        for (j = 0; j < KBYTES; j++) {
            x ^= x << 13; x ^= x >> 17; x ^= x << 5;
            m[i].kmer[j] = (uint8_t) (x & 0xff);
        }
        // force repeated 3-byte prefixes: every 5th record reuses the previous
        if (i > 0 && (i % 5) == 0) memcpy(m[i].kmer, m[i - 1].kmer, 3);
        m[i].position = (int64_t) (i * 7 % 100000);
        m[i].contig   = (int) (i % 40);
        m[i].strand   = (int) (i & 1);
    }
    qsort(m, (size_t) n, sizeof(RefRec), cmp_rec);

    rb_alloc_exact(buf, SZ, n);
    for (i = 0; i < n; i++)
        rb_set(buf, i, m[i].kmer, m[i].position, m[i].contig, m[i].strand,
               (uint8_t) (i % 41));            // distinguishable LCP bytes
    return m;
}

// K-01 -- index[] is a CUMULATIVE "prefix <= p" count, not a per-prefix count.
//
// This is the easiest thing in the adapter to get wrong, and a mistake would
// not crash: FastGA would just read the wrong slice of entries for each
// k-mer.  So the expectation is recomputed here directly from the INPUT
// records -- histogram the prefixes, then prefix-sum -- and never read back
// through the adapter.
//
// Also checked: cumulative_count_by_prefix is non-decreasing (a cumulative
// count must be), its last entry equals the total record count, and
// P.largest_prefix_group_num_entries equals the longest run of one prefix.
// It matters because the merge uses it (+1 for the end-of-group marker) to
// size a cache; too small and FastGA would overrun it.
static void k01_index_semantics(void)
{
    enum { N = 4000 };
    RecordBuffer buf;
    KmerStreamAdapter T; PayloadLayout P;
    RefRec  *m = build_sorted(&buf, N, 24601u);
    int64_t *want;
    int64_t  i, cum, run, maxrun;
    int      p, prev;

    memset(&T, 0, sizeof T); memset(&P, 0, sizeof P);
    CHECK(build_kmer_adapter(&buf, &T, &P) == 0);

    want = calloc((size_t) T.num_prefixes, sizeof(int64_t));
    for (i = 0; i < N; i++) want[prefix_of(m[i].kmer)]++;
    cum = 0;
    for (p = 0; p < T.num_prefixes; p++) { cum += want[p]; want[p] = cum; }

    CHECK_MSG(memcmp(want, T.cumulative_count_by_prefix, (size_t) T.num_prefixes * sizeof(int64_t)) == 0,
              "cumulative_count_by_prefix does not match an independent prefix<=p cumulative count");
    CHECK_EQ_I(T.cumulative_count_by_prefix[T.num_prefixes - 1], N);
    for (p = 1; p < T.num_prefixes; p++)
        if (T.cumulative_count_by_prefix[p] < T.cumulative_count_by_prefix[p - 1]) {
            t_checks++;
            t_fail_at(__FILE__, __LINE__, "cumulative_count_by_prefix non-decreasing",
                      "cumulative_count_by_prefix[%d] < [%d]", p, p - 1);
            break;
        }
    CHECK(1);

    // largest_prefix_group_num_entries: longest run of one 3-byte prefix
    maxrun = 0; run = 0; prev = -1;
    for (i = 0; i < N; i++) {
        int pr = prefix_of(m[i].kmer);
        if (pr == prev) run++; else { run = 1; prev = pr; }
        if (run > maxrun) maxrun = run;
    }
    CHECK_MSG(P.largest_prefix_group_num_entries == maxrun,
              "largest_prefix_group_num_entries = %lld, want %lld (longest run)",
              (long long) P.largest_prefix_group_num_entries, (long long) maxrun);

    free(want); free(m);
    free_kmer_adapter(&T, &P);
    rb_free(&buf);
}

// K-02 -- The derived byte widths are what the FastGA code expects.
//
// These are pure arithmetic, but every one of them is used as a stride or an
// offset by code that was written against FastK's layout, so a wrong value
// silently misaligns every entry:
//
//   prefix_bytes = 3                 prefix bytes, implied by count position
//   kmer_bytes = 10                  bytes in a 40-mer (4 bases per byte)
//   entry_kmer_suffix_bytes = kmer_bytes - prefix_bytes = 7
//                                    suffix bytes stored per entry
//   entry_bytes = suffix + 1 + 1 + pos + cont
//               = 7 + mask + lcp + 3 + 1 = 13   bytes actually stored per entry
//   prefix_plus_entry_bytes = prefix_bytes + entry_bytes = 16
//                                    width of a fully reassembled entry
//   num_prefixes = 1 << (8 * prefix_bytes) = 2^24   one count per possible prefix
static void k02_widths(void)
{
    RecordBuffer buf;
    KmerStreamAdapter T; PayloadLayout P;
    RefRec *m = build_sorted(&buf, 64, 99u);

    memset(&T, 0, sizeof T); memset(&P, 0, sizeof P);
    CHECK(build_kmer_adapter(&buf, &T, &P) == 0);

    CHECK_EQ_I(T.kmer_length_bases,  40);
    CHECK_EQ_I(T.prefix_bytes,  3);
    CHECK_EQ_I(T.kmer_bytes, 10);
    CHECK_EQ_I(T.entry_kmer_suffix_bytes,  7);           // kmer_bytes - prefix_bytes
    CHECK_EQ_I(T.entry_bytes, 7 + 1 + 1 + 3 + 1);        // suffix + mask + lcp + pos + cont
    CHECK_EQ_I(T.entry_bytes, 13);
    CHECK_EQ_I(T.prefix_plus_entry_bytes, T.prefix_bytes + T.entry_bytes);
    CHECK_EQ_I(T.prefix_plus_entry_bytes, 16);
    CHECK_EQ_I(T.num_prefixes, 1 << 24);
    CHECK_EQ_I(T.num_entries,  buf.count);
    CHECK_EQ_I(P.payload_bytes, SZ.position_bytes + SZ.contig_strand_bytes);
    CHECK_EQ_I(P.contig_strand_bytes, SZ.contig_strand_bytes);

    free(m); free_kmer_adapter(&T, &P); rb_free(&buf);
}

// K-03 -- A stored entry has exactly the layout FastGA reads.
//
// One record with distinctive values, so each field is identifiable:
//
//   bytes 0..6   the k-mer's SUFFIX (bytes 3..9); bytes 0..2 are not stored,
//                they are implied by the index position
//   byte 7       mask byte, always 0 (masking is unsupported in this pipeline)
//   byte 8       the record's LCP byte, carried through unchanged -- the merge
//                reads it to bound adaptamer groups, so dropping it would
//                silently change which seeds are found
//   bytes 9..    position and contig, copied verbatim from the record
static void k03_entry_layout(void)
{
    RecordBuffer buf;
    KmerStreamAdapter T; PayloadLayout P;
    uint8_t k[KBYTES];
    int     i;

    for (i = 0; i < KBYTES; i++) k[i] = (uint8_t) (0x10 + i);
    rb_alloc_exact(&buf, SZ, 1);
    rb_set(&buf, 0, k, 0x123456, 9, 1, 27);

    memset(&T, 0, sizeof T); memset(&P, 0, sizeof P);
    CHECK(build_kmer_adapter(&buf, &T, &P) == 0);

    // 7 suffix bytes = k-mer bytes 3..9
    CHECK_MSG(memcmp(T.entries, k + 3, 7) == 0, "suffix bytes wrong");
    CHECK_EQ_I(T.entries[7], 0);                 // mask
    CHECK_EQ_I(T.entries[8], 27);                // the record's LCP byte
    // then post+cont bytes verbatim from POS_OFFSET
    CHECK_MSG(memcmp(T.entries + 9, buf.data + POS_OFFSET,
                     (size_t) (SZ.position_bytes + SZ.contig_strand_bytes)) == 0,
              "position/contig bytes not copied verbatim");

    free_kmer_adapter(&T, &P); rb_free(&buf);
}

// K-04 -- Walking the stream visits every entry in order and reconstructs the
//          original k-mers exactly.  The strongest single adapter test.
//
// This exercises the whole addressing scheme end to end.  For each entry it
// checks that:
//   - cidx counts up one per entry, so nothing is skipped or repeated;
//   - cpre equals the TRUE prefix of that entry -- i.e. the index-walking in
//     Next_Kmer_Entry stays in step with the entries, including across gaps
//     where a prefix has no records at all;
//   - Current_Entry's first 10 bytes equal the original 40-mer, which only
//     works if the 3 implied prefix bytes and the 7 stored suffix bytes are
//     reassembled correctly.
//
// Finally it checks the end-of-stream signal (csuf == NULL and cpre == ixlen),
// since the merge loop uses that to terminate.
static void k04_traversal_roundtrip(void)
{
    enum { N = 1000 };
    RecordBuffer buf;
    KmerStreamAdapter T; PayloadLayout P;
    RefRec  *m = build_sorted(&buf, N, 555u);
    uint8_t *ent;
    int64_t  i;

    memset(&T, 0, sizeof T); memset(&P, 0, sizeof P);
    CHECK(build_kmer_adapter(&buf, &T, &P) == 0);
    ent = malloc((size_t) T.prefix_plus_entry_bytes);

    adapter_First_Kmer_Entry(&T);
    for (i = 0; i < N; i++) {
        CHECK_MSG(T.current_entry != NULL, "stream ended early at %lld", (long long) i);
        if (T.current_entry == NULL) break;
        CHECK_MSG(T.current_entry_index == i, "current_entry_index = %lld, want %lld", (long long) T.current_entry_index, (long long) i);
        CHECK_MSG(T.current_prefix == prefix_of(m[i].kmer),
                  "current_prefix = %d at record %lld, want %d",
                  T.current_prefix, (long long) i, prefix_of(m[i].kmer));
        adapter_Current_Entry(&T, ent);
        CHECK_MSG(memcmp(ent, m[i].kmer, KBYTES) == 0,
                  "Current_Entry k-mer mismatch at record %lld", (long long) i);
        adapter_Next_Kmer_Entry(&T);
    }
    // end of stream
    CHECK_MSG(T.current_entry == NULL, "stream should be exhausted after %d entries", N);
    CHECK_EQ_I(T.current_prefix, T.num_prefixes);

    free(ent); free(m); free_kmer_adapter(&T, &P); rb_free(&buf);
}

// K-05 -- Clones share the big allocations; freeing one must not free them.
//
// FastGA's merge gives each thread its own cursor over the same table by
// cloning the stream.  A clone therefore copies the small cursor state but
// SHARES the 128 MiB index and the entry table with its parent.
//
// The ownership rule that follows: adapter_Free_Kmer_Stream frees only the
// clone's wrapper, and free_kmer_adapter on the parent releases the shared
// storage.  Get that wrong and you get either a double free or a leak of
// 128 MiB per thread.
//
// Checked: the clone's table/index POINTERS are equal to the parent's (shared,
// not copied); advancing the clone does not move the parent's cursor; and
// after freeing the clone the parent is still usable -- which it would not be
// if the shared storage had been released.  Run under ASan, where a double
// free would be caught outright.
static void k05_clone_sharing(void)
{
    enum { N = 200 };
    RecordBuffer buf;
    KmerStreamAdapter T, *c; PayloadLayout P;
    RefRec *m = build_sorted(&buf, N, 31u);
    int64_t parent_cidx;

    memset(&T, 0, sizeof T); memset(&P, 0, sizeof P);
    CHECK(build_kmer_adapter(&buf, &T, &P) == 0);

    adapter_First_Kmer_Entry(&T);
    parent_cidx = T.current_entry_index;

    c = adapter_Clone_Kmer_Stream(&T);
    CHECK_MSG(c->entries == T.entries, "clone does not share entries");
    CHECK_MSG(c->cumulative_count_by_prefix == T.cumulative_count_by_prefix,
              "clone does not share cumulative_count_by_prefix");
    CHECK_EQ_I(c->is_clone, 1);
    CHECK_EQ_I(T.is_clone, 0);

    adapter_Next_Kmer_Entry(c);
    adapter_Next_Kmer_Entry(c);
    CHECK_MSG(T.current_entry_index == parent_cidx, "advancing the clone moved the parent");
    CHECK_MSG(c->current_entry_index == parent_cidx + 2, "clone did not advance");

    adapter_Free_Kmer_Stream(c);        // frees the wrapper only
    // parent still usable: shared storage must not have been released
    adapter_Next_Kmer_Entry(&T);
    CHECK_MSG(T.current_entry_index == parent_cidx + 1, "parent unusable after freeing a clone");

    free(m); free_kmer_adapter(&T, &P); rb_free(&buf);
}

int main(void)
{
    SZ = rs_make(3, 1);
    printf("test_adapter\n");
    RUN(k01_index_semantics);
    RUN(k02_widths);
    RUN(k03_entry_layout);
    RUN(k04_traversal_roundtrip);
    RUN(k05_clone_sharing);
    return t_report("test_adapter");
}
