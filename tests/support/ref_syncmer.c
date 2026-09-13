#include <stdlib.h>
#include <string.h>
#include "ref_syncmer.h"
#include "synth_gdb.h"
#include "tables.h"      // TMap[] only: constant data, part of the definition

#define REF_KMER 40
#define REF_TMER 12
#define REF_SMER  8
#define REF_SOFF  (REF_TMER - REF_SMER)

void ref_list_init(RefRecList *l) { memset(l, 0, sizeof(*l)); }
void ref_list_free(RefRecList *l) { free(l->v); memset(l, 0, sizeof(*l)); }

static void ref_push(RefRecList *l, int contig, int64_t pos, int strand,
                     const uint8_t *kmer)
{
    if (l->n == l->cap) {
        l->cap = l->cap ? l->cap * 2 : 1024;
        l->v   = (RefRec *) realloc(l->v, (size_t) l->cap * sizeof(RefRec));
    }
    l->v[l->n].contig   = contig;
    l->v[l->n].position = pos;
    l->v[l->n].strand   = strand;
    memcpy(l->v[l->n].kmer, kmer, 10);
    l->n++;
}

// Pack 4 bases, MSB first.
static int ref_pack4(const char *seq, int64_t at)
{
    return (base_code(seq[at])     << 6) | (base_code(seq[at + 1]) << 4)
         | (base_code(seq[at + 2]) << 2) |  base_code(seq[at + 3]);
}

// Reverse complement of a packed 4-base byte, from first principles.
//
// Two things happen at once: the four bases REVERSE order, and each is
// COMPLEMENTED.  With the encoding A=0 C=1 G=2 T=3, complementing is just
// 3 - v (A<->T is 0<->3, C<->G is 1<->2).  Field i counting from the low end
// moves to position 3-i counting from the high end.
//
// Production does this with the precomputed Comp[] table; recomputing it here
// keeps the oracle genuinely independent.
static int ref_rc4(int b)
{
    int out = 0, i;
    for (i = 0; i < 4; i++) {
        int v = (b >> (2 * i)) & 3;          // field i counting from the low end
        out |= (3 - v) << (2 * (3 - i));
    }
    return out & 0xff;
}

// Canonical hash of the 8-mer at offset i.
//
// The 8-mer spans two packed bytes: bases i..i+3 (lo) and i+4..i+7 (hi).  Its
// forward hash scrambles each byte through TMap and concatenates them, high
// half first.  Its reverse-complement hash does the same to the RC, where the
// two bytes both swap places and get complemented.
//
// Taking the MINIMUM of the two makes the value strand-independent, so the
// same physical position is selected whichever strand you happen to read.
static int ref_mz(const char *seq, int64_t i)
{
    int lo  = ref_pack4(seq, i);
    int hi  = ref_pack4(seq, i + 4);
    int fwd = (TMap[lo] << 8) | TMap[hi];
    int rc  = (TMap[ref_rc4(hi)] << 8) | TMap[ref_rc4(lo)];
    return (fwd < rc) ? fwd : rc;
}

void ref_pack_forward(const char *seq, int64_t j, uint8_t *out)
{
    int b;
    for (b = 0; b < 10; b++) out[b] = (uint8_t) ref_pack4(seq, j + 4 * b);
}

// Pack the reverse-complement 40-mer keyed at offset j.
//
// It covers FORWARD bases [j-28, j+12), reverse-complemented.  The asymmetric
// window is why the RC record's stored position is j+12 rather than j: the
// k-mer's own start corresponds to the right-hand end of that forward span.
//
// Built explicitly as a complemented, reversed base string and then packed,
// rather than by manipulating packed bytes, so it is obviously correct by
// inspection.
void ref_pack_rc(const char *seq, int64_t j, uint8_t *out)
{
    // RC of forward [j-28, j+12): RC[i] = complement(seq[j + 11 - i]).
    char rc[REF_KMER + 1];
    int  i, b;
    for (i = 0; i < REF_KMER; i++)
        rc[i] = base_char(3 - base_code(seq[j + 11 - i]));
    rc[REF_KMER] = '\0';
    for (b = 0; b < 10; b++) out[b] = (uint8_t) ref_pack4(rc, 4 * b);
}

int ref_is_selected(const char *seq, int len, int64_t j)
{
    int h[5], m, a, t;
    if (j < 0 || j > len - REF_TMER) return 0;
    for (t = 0; t < 5; t++) h[t] = ref_mz(seq, j + t);
    m = h[0]; a = 0;
    for (t = 1; t < 5; t++) if (h[t] < m) { m = h[t]; a = t; }
    return (a == 0 || h[4] == m);
}

void ref_extract_contig(const char *seq, int len, int contig, RefRecList *out)
{
    int64_t j;
    if (len < REF_TMER) return;

    for (j = 0; j <= len - REF_TMER; j++) {
        int h[5], m, a, t;

        for (t = 0; t < 5; t++) h[t] = ref_mz(seq, j + t);

        m = h[0]; a = 0;
        for (t = 1; t < 5; t++) if (h[t] < m) { m = h[t]; a = t; }   // earliest argmin

        if (a == 0)                      out->n_left++;
        if (h[4] < m || (h[4] == m && a == 4)) out->n_right++;
        if (h[4] == m && a < 4)          out->n_tie++;

        if (!(a == 0 || h[4] == m)) continue;                        // not a closed syncmer
        out->n_sel++;

        if (j <= len - REF_KMER) {
            uint8_t k[10];
            ref_pack_forward(seq, j, k);
            ref_push(out, contig, j, 0, k);
        }
        if (j >= REF_KMER - REF_TMER) {
            uint8_t k[10];
            ref_pack_rc(seq, j, k);
            ref_push(out, contig, j + REF_TMER, 1, k);
        }
    }
}

// Independent base-level LCP: how many leading BASES two packed k-mers share.
//
// Unpacks both to plain base arrays and counts the common prefix with a simple
// loop.  Slow and obvious on purpose -- production's simple_kmer_lcp_bases
// works on packed bytes and has to locate the differing base inside a byte,
// which is exactly the kind of bit manipulation worth checking against
// something dumber.  Used anywhere an LCP expectation is needed, including the
// MPI seam tests, so those checks are never self-referential.
int ref_lcp_bases(const uint8_t *a, const uint8_t *b)
{
    char sa[REF_KMER], sb[REF_KMER];
    int  i, bb;
    for (bb = 0; bb < 10; bb++) {
        for (i = 0; i < 4; i++) {
            sa[bb * 4 + i] = (char) ((a[bb] >> (2 * (3 - i))) & 3);
            sb[bb * 4 + i] = (char) ((b[bb] >> (2 * (3 - i))) & 3);
        }
    }
    for (i = 0; i < REF_KMER; i++) if (sa[i] != sb[i]) return i;
    return REF_KMER;
}

// Total order over records, used to compare two sets as MULTISETS: sort both
// and walk them in step.  Ordering by k-mer first and then by the full payload
// means equal k-mers still compare deterministically, so the comparison does
// not depend on the order records happened to arrive in.
int ref_rec_cmp(const void *x, const void *y)
{
    const RefRec *p = (const RefRec *) x, *q = (const RefRec *) y;
    int c = memcmp(p->kmer, q->kmer, 10);
    if (c) return c;
    if (p->position != q->position) return (p->position < q->position) ? -1 : 1;
    if (p->contig   != q->contig)   return (p->contig   < q->contig)   ? -1 : 1;
    return p->strand - q->strand;
}

void ref_list_sort(RefRecList *l)
{
    qsort(l->v, (size_t) l->n, sizeof(RefRec), ref_rec_cmp);
}
