// Read GIXshow output on stdin, emit "<contig> <position> <strand> <lcp> <hex_kmer>\n" on stdout.
//
// GIXshow LCP encoding (per inspection of yeast7.gix):
//   "0"   for the very first record (LCP undefined → starts at 0)
//   "*"   for records identical to predecessor (LCP = KMER = 40)
//   N     for byte-resolution LCP values 1..KMER-1 (= "shared bases")
//
// We normalize to a single integer in the LCP column. Sort key for the
// canonical comparison is then (contig, position, strand, kmer, lcp), but
// since LCP is fully determined by the sorted-table position rather than by
// the record itself, we keep it as a separate sortable column to preserve the
// stable ordering verify scripts use.
//
// IMPORTANT: GIXshow's LCP column is computed assuming records are walked in
// the .gix file's sort order. If you sort the normalized output by anything
// other than that order, the LCP values stop matching what msd_sort produced.
// Use this gold artifact only for cmp against pga-mpi output that has been
// sorted by the same canonical key, or compare in raw GIXshow order.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define KMER 40

static const unsigned char base_to_bits[256] = {
    ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
    ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3,
};

static const char hex_chars[] = "0123456789abcdef";

int main(void) {
    char line[512];
    char out[256];
    while (fgets(line, sizeof(line), stdin)) {
        char *colon = strchr(line, ':');
        if (!colon) continue;
        char *p = colon + 1;
        while (*p == ' ') p++;
        char *kmer = p;
        if (kmer[0] != 'a' && kmer[0] != 'c' && kmer[0] != 'g' && kmer[0] != 't' &&
            kmer[0] != 'A' && kmer[0] != 'C' && kmer[0] != 'G' && kmer[0] != 'T') continue;

        // Pack the k-mer into 10 bytes.
        unsigned char kbytes[10];
        for (int i = 0; i < 10; i++) {
            unsigned int b = 0;
            for (int k = 0; k < 4; k++) {
                b = (b << 2) | base_to_bits[(unsigned char)kmer[i*4+k]];
            }
            kbytes[i] = (unsigned char)b;
        }

        // Tokenize the rest: skip kmer, skip mask, read lcp, read sign, contig, '|', position.
        char *q = kmer + KMER;
        // mask
        while (*q == ' ') q++;
        while (*q && *q != ' ' && *q != '\n') q++;
        // lcp
        while (*q == ' ') q++;
        char *lcp_start = q;
        while (*q && *q != ' ' && *q != '\n') q++;
        char lcp_save = *q;
        if (*q) *q++ = '\0';
        // sign
        while (*q == ' ') q++;
        char sign = *q;
        while (*q && *q != ' ' && *q != '\n') q++;
        // contig
        while (*q == ' ') q++;
        char *contig_start = q;
        while (*q && *q != ' ' && *q != '\n') q++;
        if (*q) *q++ = '\0';
        // skip '|'
        while (*q == ' ') q++;
        if (*q == '|') q++;
        // position
        while (*q == ' ') q++;
        char *pos_start = q;
        while (*q && *q != ' ' && *q != '\n') q++;
        if (*q) *q = '\0';

        (void) lcp_save;

        // Decode LCP: "*" -> KMER, "0" -> 0, otherwise integer parse.
        int lcp;
        if (lcp_start[0] == '*' && lcp_start[1] == '\0')
            lcp = KMER;
        else
            lcp = atoi(lcp_start);

        int strand = (sign == '+') ? 0 : 1;

        int n = snprintf(out, sizeof(out), "%s %s %d %d ",
                         contig_start, pos_start, strand, lcp);
        for (int i = 0; i < 10; i++) {
            out[n++] = hex_chars[kbytes[i] >> 4];
            out[n++] = hex_chars[kbytes[i] & 0xf];
        }
        out[n++] = '\n';
        fwrite(out, 1, n, stdout);
    }
    return 0;
}
