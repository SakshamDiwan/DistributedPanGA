#include <stdlib.h>
#include <string.h>
#include "synth_gdb.h"
#include "gene_core.h"

int base_code(char c)
{
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return -1;
    }
}

char base_char(int v)
{
    static const char s[4] = { 'A', 'C', 'G', 'T' };
    return (v >= 0 && v < 4) ? s[v] : '?';
}

// A deterministic pseudorandom ACGT string.  xorshift32 rather than rand() so
// fixtures are reproducible across machines and libc versions -- when a test
// fails, the seed printed in its output regenerates the exact input.
char *synth_random_seq(int len, uint32_t seed)
{
    char    *s = (char *) malloc((size_t) len + 1);
    uint32_t x = seed ? seed : 1u;
    int      i;

    if (s == NULL) return NULL;
    for (i = 0; i < len; i++) {
        x ^= x << 13; x ^= x >> 17; x ^= x << 5;   // xorshift32
        s[i] = base_char((int) (x & 3u));
    }
    s[len] = '\0';
    return s;
}

int synth_gdb_build(SynthGDB *sg, const char *const *seqs, int nseq)
{
    int64_t total = 0, maxc = 0;
    size_t  need  = 0;
    int     i;

    memset(sg, 0, sizeof(*sg));
    if (nseq <= 0) return -1;

    for (i = 0; i < nseq; i++) {
        int64_t n = (int64_t) strlen(seqs[i]);
        total += n;
        if (n > maxc) maxc = n;
        need += (size_t) COMPRESSED_LEN(n);
    }

    sg->bps_len = need ? need : 1;
    sg->bps     = (uint8_t *) calloc(sg->bps_len, 1);
    sg->gdb.contigs = (GDB_CONTIG *) calloc((size_t) nseq + 1, sizeof(GDB_CONTIG));
    if (sg->bps == NULL || sg->gdb.contigs == NULL) { synth_gdb_free(sg); return -1; }

    // Lay the contigs out exactly as a .bps file does: each one starts at a
    // BYTE-ALIGNED offset and occupies COMPRESSED_LEN(clen) = (clen+3)/4
    // bytes.  There is no bit packing across a contig boundary, so boff is
    // simply the running byte total.
    //
    // The 2-bit packing itself is done with the library's own Compress_Read
    // rather than a hand-rolled packer.  That is deliberate: this is a FILE
    // FORMAT, not the logic under test, and using the library's compressor
    // guarantees its Uncompress_Read inverts it exactly.  The independence
    // that matters for these tests is in ref_syncmer.c, which re-derives
    // which k-mers should be emitted.
    {
        int64_t off = 0;
        for (i = 0; i < nseq; i++) {
            int64_t n    = (int64_t) strlen(seqs[i]);
            int     nby  = COMPRESSED_LEN(n);
            char   *tmp  = (char *) malloc((size_t) n + 4);
            int64_t k;

            if (tmp == NULL) { synth_gdb_free(sg); return -1; }
            for (k = 0; k < n; k++) {
                int c = base_code(seqs[i][k]);
                if (c < 0) { free(tmp); synth_gdb_free(sg); return -1; }
                tmp[k] = (char) c;
            }
            Compress_Read((int) n, tmp);
            memcpy(sg->bps + off, tmp, (size_t) nby);
            free(tmp);

            sg->gdb.contigs[i].clen = n;
            sg->gdb.contigs[i].boff = off;
            sg->gdb.contigs[i].sbeg = 0;
            sg->gdb.contigs[i].scaf = 0;
            off += nby;
        }
        sg->gdb.contigs[nseq].boff = off;   // terminal entry, as Read_GDB writes
    }

    sg->gdb.ncontig  = nseq;
    sg->gdb.maxctg   = maxc;
    sg->gdb.seqtot   = total;
    sg->gdb.seqstate = EXTERNAL;
    sg->gdb.seqs     = fmemopen(sg->bps, sg->bps_len, "rb");
    if (sg->gdb.seqs == NULL) { synth_gdb_free(sg); return -1; }
    return 0;
}

void synth_gdb_free(SynthGDB *sg)
{
    // Deliberately NOT Close_GDB.  That function frees headers, scaffolds,
    // srcpath and seqpath, none of which this fixture ever allocated, and
    // would fclose/free seqs under rules that depend on seqstate.  We free
    // exactly what we own: the stream and the two arrays.
    if (sg->gdb.seqs != NULL) { fclose((FILE *) sg->gdb.seqs); sg->gdb.seqs = NULL; }
    free(sg->gdb.contigs); sg->gdb.contigs = NULL;
    free(sg->bps);         sg->bps = NULL;
    sg->bps_len = 0;
}
