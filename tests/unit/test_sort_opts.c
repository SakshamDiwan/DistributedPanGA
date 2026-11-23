// ============================================================================
// R3 -- sort_records thread count and part granularity.  Tests SO-01..SO-05.
// ============================================================================
//
// sort_records_opts(buf, threads, parts) regroups records by part (first k-mer
// byte for 256 parts, first five bases for 1024) and hands the parts to
// msd_sort, one thread per part.  Changing either setting must not change
// what comes out:
//
//   * the records are sorted by k-mer and none is lost or duplicated;
//   * after recalc_all_lcps, the (k-mer, LCP) sequence is identical to the
//     default (4 threads, 256 parts) -- that sequence is all the downstream
//     merge sees apart from the payload of equal k-mers, whose relative order
//     msd_sort never specified.
//
// Fixtures are skewed on purpose: most records share a few leading bytes (so
// one 256-part holds most of the work and the 1024 split matters), there are
// many duplicate k-mers, and the extremes 0x00.. and 0xff.. are present.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "sort.h"
#include "extract.h"
#include "record.h"
#include "constants.h"

static RecordSizing SZ = { 3, 1, 16 };

static uint64_t rng = 88172645463325252ULL;
static uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

// n records; `skew` in [0,100]: percent drawn from 3 hot leading bytes.
static void fill(RecordBuffer *b, int64_t n, int skew)
{
    record_buffer_init(b, SZ);
    while (b->capacity < n + 1) {   // grow without touching production growth
        free(b->data);
        b->capacity *= 2;
        b->data = malloc(b->capacity * SZ.record_size + 1);
    }
    static const uint8_t hot[3] = { 0x00, 0x3f, 0xff };
    for (int64_t i = 0; i < n; i++) {
        uint8_t *r = b->data + i * SZ.record_size;
        uint64_t v = xr();
        r[0] = (uint8_t) v;                          // garbage LCP byte
        for (int k = 0; k < KBYTES; k++) r[1 + k] = (uint8_t) xr();
        if ((int) (v % 100) < skew) {
            r[1] = hot[(v >> 8) % 3];
            r[2] &= (v >> 12) & 1 ? 0x3f : 0xff;     // also skew the 5th base
            if ((v >> 16) % 4 == 0) memset(r + 4, 0, KBYTES - 3);   // long LCPs
        }
        if (i > 0 && (v >> 20) % 5 == 0)             // duplicate k-mer, new payload
            memcpy(r + 1, b->data + (i - 1 - (int64_t) ((v >> 24) % (uint64_t) (i < 8 ? i : 8))) * SZ.record_size + 1, KBYTES);
        r[MASK_OFFSET] = 0;
        for (int k = POS_OFFSET; k < SZ.record_size; k++) r[k] = (uint8_t) xr();
    }
    b->count = n;
}

static RecordBuffer clone(const RecordBuffer *a)
{
    RecordBuffer b = *a;
    b.data = malloc(a->count * SZ.record_size + 1);
    memcpy(b.data, a->data, a->count * SZ.record_size);
    b.capacity = a->count;
    return b;
}

static int cmp_full(const void *x, const void *y)   // whole record except LCP byte
{
    return memcmp((const uint8_t *) x + 1, (const uint8_t *) y + 1, SZ.record_size - 1);
}

// Sorted (by full record) copy, for multiset comparison.
static uint8_t *multiset(const RecordBuffer *b)
{
    uint8_t *m = malloc(b->count * SZ.record_size + 1);
    memcpy(m, b->data, b->count * SZ.record_size);
    qsort(m, b->count, SZ.record_size, cmp_full);
    return m;
}

static int same_multiset(const RecordBuffer *a, const RecordBuffer *b)
{
    if (a->count != b->count) return 0;
    uint8_t *x = multiset(a), *y = multiset(b);
    int ok = 1;
    for (int64_t i = 0; i < a->count && ok; i++)
        ok = memcmp(x + i * SZ.record_size + 1, y + i * SZ.record_size + 1, SZ.record_size - 1) == 0;
    free(x); free(y);
    return ok;
}

static int same_kmer_lcp_sequence(const RecordBuffer *a, const RecordBuffer *b)
{
    if (a->count != b->count) return 0;
    for (int64_t i = 0; i < a->count; i++)
        if (memcmp(a->data + i * SZ.record_size, b->data + i * SZ.record_size, 1 + KBYTES) != 0)
            return 0;
    return 1;
}

static int is_sorted(const RecordBuffer *b)
{
    for (int64_t i = 1; i < b->count; i++)
        if (memcmp(b->data + (i - 1) * SZ.record_size + 1, b->data + i * SZ.record_size + 1, KBYTES) > 0)
            return 0;
    return 1;
}

// One configuration against the default; returns through CHECKs.
static void against_default(const RecordBuffer *in, int threads, int parts, const char *what)
{
    RecordBuffer ref = clone(in), got = clone(in);
    sort_records_opts(&ref, 4, 256);
    recalc_all_lcps(&ref);
    sort_records_opts(&got, threads, parts);
    CHECK_MSG(g_sort_timing.threads == threads && g_sort_timing.parts == parts,
              "%s: timing reports threads=%d parts=%d", what, g_sort_timing.threads, g_sort_timing.parts);
    recalc_all_lcps(&got);
    CHECK_MSG(is_sorted(&got), "%s: threads=%d parts=%d not sorted", what, threads, parts);
    CHECK_MSG(same_multiset(in, &got), "%s: threads=%d parts=%d changed the record multiset", what, threads, parts);
    CHECK_MSG(same_kmer_lcp_sequence(&ref, &got), "%s: threads=%d parts=%d (k-mer, LCP) sequence differs from default",
              what, threads, parts);
    free(ref.data); free(got.data);
}

// SO-01 -- every thread count x part setting gives the default's (k-mer, LCP)
// sequence and the input multiset, on uniform and on heavily skewed input.
static void so01_settings_equivalent(void)
{
    static const int T[] = { 1, 2, 4, 7, 16 };
    for (int skew = 0; skew <= 95; skew += 95) {
        RecordBuffer in;
        fill(&in, 60000, skew);
        for (int t = 0; t < 5; t++)
            for (int p = 256; p <= 1024; p += 768)
                against_default(&in, T[t], p, skew ? "skewed" : "uniform");
        record_buffer_free(&in);
    }
}

// SO-02 -- degenerate sizes: 0, 1, 2 records, and everything in one 5-base
// bucket (one part with 1024, one part with 256 too).
static void so02_degenerate(void)
{
    for (int n = 0; n <= 2; n++) {
        RecordBuffer in;
        fill(&in, n, 50);
        against_default(&in, 4, 1024, "tiny");
        record_buffer_free(&in);
    }
    RecordBuffer in;
    fill(&in, 20000, 0);
    for (int64_t i = 0; i < in.count; i++) {     // all in bucket 1023 (TTTTT...)
        in.data[i * SZ.record_size + 1] = 0xff;
        in.data[i * SZ.record_size + 2] |= 0xc0;
    }
    against_default(&in, 4, 1024, "one bucket");
    CHECK_EQ_I(g_sort_timing.nonempty_parts, 1);
    CHECK_EQ_I(g_sort_timing.largest_part, 20000);
    record_buffer_free(&in);
}

// SO-03 -- nonempty_parts counts the distinct parts present.
static void so03_part_counts(void)
{
    RecordBuffer in;
    fill(&in, 4, 0);
    const uint8_t first[4][2] = { { 0x00, 0x00 }, { 0x00, 0x40 }, { 0x00, 0x41 }, { 0x01, 0x00 } };
    for (int i = 0; i < 4; i++) memcpy(in.data + i * SZ.record_size + 1, first[i], 2);
    RecordBuffer a = clone(&in), b = clone(&in);
    sort_records_opts(&a, 4, 256);
    CHECK_EQ_I(g_sort_timing.nonempty_parts, 2);      // bytes 0x00, 0x01
    sort_records_opts(&b, 4, 1024);
    CHECK_EQ_I(g_sort_timing.nonempty_parts, 3);      // buckets 0, 1, 4
    free(a.data); free(b.data); record_buffer_free(&in);
}

// SO-04 -- environment parsing and defaults.
static void so04_env(void)
{
    unsetenv("PGA_SORT_THREADS"); unsetenv("PGA_SORT_PARTS");
    CHECK_EQ_I(sort_threads_from_env(), 4);
    CHECK_EQ_I(sort_parts_from_env(), 256);
    setenv("PGA_SORT_THREADS", "8", 1);   CHECK_EQ_I(sort_threads_from_env(), 8);
    setenv("PGA_SORT_THREADS", "0", 1);   CHECK_EQ_I(sort_threads_from_env(), 4);
    setenv("PGA_SORT_THREADS", "999", 1); CHECK_EQ_I(sort_threads_from_env(), 4);
    setenv("PGA_SORT_THREADS", "x", 1);   CHECK_EQ_I(sort_threads_from_env(), 4);
    setenv("PGA_SORT_PARTS", "1024", 1);  CHECK_EQ_I(sort_parts_from_env(), 1024);
    setenv("PGA_SORT_PARTS", "4096", 1);  CHECK_EQ_I(sort_parts_from_env(), 256);
    setenv("PGA_SORT_PARTS", "", 1);      CHECK_EQ_I(sort_parts_from_env(), 256);
    unsetenv("PGA_SORT_THREADS"); unsetenv("PGA_SORT_PARTS");
}

// SO-05 -- sort_records() with the environment unset is the 4/256 default,
// and with it set follows it.
static void so05_sort_records_uses_env(void)
{
    RecordBuffer in;
    fill(&in, 5000, 50);
    RecordBuffer a = clone(&in);
    unsetenv("PGA_SORT_THREADS"); unsetenv("PGA_SORT_PARTS");
    sort_records(&a);
    CHECK_EQ_I(g_sort_timing.threads, 4);
    CHECK_EQ_I(g_sort_timing.parts, 256);
    RecordBuffer b = clone(&in);
    setenv("PGA_SORT_THREADS", "3", 1); setenv("PGA_SORT_PARTS", "1024", 1);
    sort_records(&b);
    CHECK_EQ_I(g_sort_timing.threads, 3);
    CHECK_EQ_I(g_sort_timing.parts, 1024);
    unsetenv("PGA_SORT_THREADS"); unsetenv("PGA_SORT_PARTS");
    recalc_all_lcps(&a); recalc_all_lcps(&b);
    CHECK(same_kmer_lcp_sequence(&a, &b));
    free(a.data); free(b.data); record_buffer_free(&in);
}

int main(void)
{
    RUN(so01_settings_equivalent);
    RUN(so02_degenerate);
    RUN(so03_part_counts);
    RUN(so04_env);
    RUN(so05_sort_records_uses_env);
    return t_report("test_sort_opts");
}
