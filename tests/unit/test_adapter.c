// ============================================================================
// P4 -- The FastK/FastGA adapter.  Tests K-01..K-05.  See tests/CATALOG.md.
// ============================================================================
//
// WHAT THE ADAPTER IS FOR.  The downstream merge and alignment code is lifted
// from FastGA and expects to read k-mers through FastK's Kmer_Stream interface.
// Our pipeline instead produces a flat array of our own records.  Rather than
// rewrite FastGA, src/kmer_adapter.c builds a structure whose FIELD LAYOUT
// MIRRORS Kmer_Stream, and the call sites simply cast to it.
//
// That cast is why these tests matter: nothing in the compiler checks the
// mirror.  Reorder a field, change a width, or get the index semantics subtly
// wrong, and FastGA reads plausible garbage instead of failing loudly.
//
// HOW THE STREAM IS ADDRESSED.  A 40-mer is 10 bytes.  The adapter splits it:
//   - the first ibyte = 3 bytes are the PREFIX.  These are not stored per
//     entry at all; the position in the index array implies them.
//   - the remaining hbyte = 7 bytes are stored with each entry, along with a
//     mask byte, the LCP byte, and the position/contig payload.
// index[p] is CUMULATIVE: the number of records whose prefix is <= p.  So the
// entries for prefix p occupy [index[p-1], index[p]), and a reader tracks the
// current prefix in `cpre` while walking entries.  Current_Entry() rebuilds
// the full k-mer by prepending cpre's three bytes to the stored suffix.
//
// SIZE WARNING: ixlen = 2^24 entries of int64 = 128 MiB per adapter, allocated
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
// RUNS of equal prefixes.  Without them index[] would be nearly all 0s and 1s,
// maxp would be trivially 1, and K-01 would not distinguish a cumulative count
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
// Also checked: index[] is non-decreasing (a cumulative count must be), its
// last entry equals the total record count, and P.maxp equals the longest run
// of one prefix plus 1.  maxp matters because the merge uses it to size a
// cache; too small and FastGA would overrun it.
static void k01_index_semantics(void)
{
    enum { N = 4000 };
    RecordBuffer buf;
    KmerStreamAdapter T; PostListAdapter P;
    RefRec  *m = build_sorted(&buf, N, 24601u);
    int64_t *want;
    int64_t  i, cum, run, maxrun;
    int      p, prev;

    memset(&T, 0, sizeof T); memset(&P, 0, sizeof P);
    CHECK(build_kmer_adapter(&buf, &T, &P) == 0);

    want = calloc((size_t) T.ixlen, sizeof(int64_t));
    for (i = 0; i < N; i++) want[prefix_of(m[i].kmer)]++;
    cum = 0;
    for (p = 0; p < T.ixlen; p++) { cum += want[p]; want[p] = cum; }

    CHECK_MSG(memcmp(want, T.index, (size_t) T.ixlen * sizeof(int64_t)) == 0,
              "index does not match an independent prefix<=p cumulative count");
    CHECK_EQ_I(T.index[T.ixlen - 1], N);
    for (p = 1; p < T.ixlen; p++)
        if (T.index[p] < T.index[p - 1]) {
            t_checks++;
            t_fail_at(__FILE__, __LINE__, "index non-decreasing",
                      "index[%d] < index[%d]", p, p - 1);
            break;
        }
    CHECK(1);

    // maxp is the longest run of one 3-byte prefix, plus 1
    maxrun = 0; run = 0; prev = -1;
    for (i = 0; i < N; i++) {
        int pr = prefix_of(m[i].kmer);
        if (pr == prev) run++; else { run = 1; prev = pr; }
        if (run > maxrun) maxrun = run;
    }
    CHECK_MSG(P.maxp == maxrun + 1, "maxp = %lld, want %lld (longest run %lld + 1)",
              (long long) P.maxp, (long long) (maxrun + 1), (long long) maxrun);

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
//   ibyte = 3                    prefix bytes, implied by index position
//   kbyte = 10                   bytes in a 40-mer (4 bases per byte)
//   hbyte = kbyte - ibyte = 7    suffix bytes stored per entry
//   pbyte = hbyte + 1 + 1 + post + cont
//         = 7 + mask + lcp + 3 + 1 = 13   bytes actually stored per entry
//   tbyte = ibyte + pbyte = 16   width of a fully reassembled entry
//   ixlen = 1 << (8 * ibyte) = 2^24   one index slot per possible prefix
static void k02_widths(void)
{
    RecordBuffer buf;
    KmerStreamAdapter T; PostListAdapter P;
    RefRec *m = build_sorted(&buf, 64, 99u);

    memset(&T, 0, sizeof T); memset(&P, 0, sizeof P);
    CHECK(build_kmer_adapter(&buf, &T, &P) == 0);

    CHECK_EQ_I(T.kmer,  40);
    CHECK_EQ_I(T.ibyte,  3);
    CHECK_EQ_I(T.kbyte, 10);
    CHECK_EQ_I(T.hbyte,  7);                       // kbyte - ibyte
    CHECK_EQ_I(T.pbyte, 7 + 1 + 1 + 3 + 1);        // hbyte + mask + lcp + pos + cont
    CHECK_EQ_I(T.pbyte, 13);
    CHECK_EQ_I(T.tbyte, T.ibyte + T.pbyte);
    CHECK_EQ_I(T.tbyte, 16);
    CHECK_EQ_I(T.ixlen, 1 << 24);
    CHECK_EQ_I(T.nels,  buf.count);
    CHECK_EQ_I(P.pbyte, SZ.post_bytes + SZ.cont_bytes);
    CHECK_EQ_I(P.cbyte, SZ.cont_bytes);

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
    KmerStreamAdapter T; PostListAdapter P;
    uint8_t k[KBYTES];
    int     i;

    for (i = 0; i < KBYTES; i++) k[i] = (uint8_t) (0x10 + i);
    rb_alloc_exact(&buf, SZ, 1);
    rb_set(&buf, 0, k, 0x123456, 9, 1, 27);

    memset(&T, 0, sizeof T); memset(&P, 0, sizeof P);
    CHECK(build_kmer_adapter(&buf, &T, &P) == 0);

    // 7 suffix bytes = k-mer bytes 3..9
    CHECK_MSG(memcmp(T.table, k + 3, 7) == 0, "suffix bytes wrong");
    CHECK_EQ_I(T.table[7], 0);                 // mask
    CHECK_EQ_I(T.table[8], 27);                // the record's LCP byte
    // then post+cont bytes verbatim from POS_OFFSET
    CHECK_MSG(memcmp(T.table + 9, buf.data + POS_OFFSET,
                     (size_t) (SZ.post_bytes + SZ.cont_bytes)) == 0,
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
    KmerStreamAdapter T; PostListAdapter P;
    RefRec  *m = build_sorted(&buf, N, 555u);
    uint8_t *ent;
    int64_t  i;

    memset(&T, 0, sizeof T); memset(&P, 0, sizeof P);
    CHECK(build_kmer_adapter(&buf, &T, &P) == 0);
    ent = malloc((size_t) T.tbyte);

    adapter_First_Kmer_Entry(&T);
    for (i = 0; i < N; i++) {
        CHECK_MSG(T.csuf != NULL, "stream ended early at %lld", (long long) i);
        if (T.csuf == NULL) break;
        CHECK_MSG(T.cidx == i, "cidx = %lld, want %lld", (long long) T.cidx, (long long) i);
        CHECK_MSG(T.cpre == prefix_of(m[i].kmer),
                  "cpre = %d at record %lld, want %d",
                  T.cpre, (long long) i, prefix_of(m[i].kmer));
        adapter_Current_Entry(&T, ent);
        CHECK_MSG(memcmp(ent, m[i].kmer, KBYTES) == 0,
                  "Current_Entry k-mer mismatch at record %lld", (long long) i);
        adapter_Next_Kmer_Entry(&T);
    }
    // end of stream
    CHECK_MSG(T.csuf == NULL, "stream should be exhausted after %d entries", N);
    CHECK_EQ_I(T.cpre, T.ixlen);

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
    KmerStreamAdapter T, *c; PostListAdapter P;
    RefRec *m = build_sorted(&buf, N, 31u);
    int64_t parent_cidx;

    memset(&T, 0, sizeof T); memset(&P, 0, sizeof P);
    CHECK(build_kmer_adapter(&buf, &T, &P) == 0);

    adapter_First_Kmer_Entry(&T);
    parent_cidx = T.cidx;

    c = adapter_Clone_Kmer_Stream(&T);
    CHECK_MSG(c->table == T.table, "clone does not share table");
    CHECK_MSG(c->index == T.index, "clone does not share index");
    CHECK_EQ_I(c->clone, 1);
    CHECK_EQ_I(T.clone, 0);

    adapter_Next_Kmer_Entry(c);
    adapter_Next_Kmer_Entry(c);
    CHECK_MSG(T.cidx == parent_cidx, "advancing the clone moved the parent");
    CHECK_MSG(c->cidx == parent_cidx + 2, "clone did not advance");

    adapter_Free_Kmer_Stream(c);        // frees the wrapper only
    // parent still usable: shared storage must not have been released
    adapter_Next_Kmer_Entry(&T);
    CHECK_MSG(T.cidx == parent_cidx + 1, "parent unusable after freeing a clone");

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
