// Read GIXshow output on stdin, emit "<contig> <position> <strand> <hex_kmer>\n" on stdout.
// GIXshow line: "<idx>: <40-base-kmer> mask lcp sign contig | position"
//   - mask is '*' or hex
//   - lcp is '*' or int
//   - sign is '+' or '-'
//   - the '|' is a literal pipe character
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

static const unsigned char base_to_bits[256] = {
    ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
    ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3,
};

static const char hex_chars[] = "0123456789abcdef";

int main(void) {
    char line[512];
    char out[256];
    while (fgets(line, sizeof(line), stdin)) {
        // Find the colon ending the index. Lines without that pattern (like the header) are skipped.
        char *colon = strchr(line, ':');
        if (!colon) continue;
        char *p = colon + 1;
        while (*p == ' ') p++;
        // Now p points at the start of the k-mer.
        char *kmer = p;
        // K-mer is exactly 40 chars (a/c/g/t).
        if (strlen(kmer) < 40 || !base_to_bits[(unsigned char)kmer[0]] && kmer[0] != 'a' && kmer[0] != 'A') {
            // Validate first base; skip if not a DNA base.
            if (kmer[0] != 'a' && kmer[0] != 'c' && kmer[0] != 'g' && kmer[0] != 't' &&
                kmer[0] != 'A' && kmer[0] != 'C' && kmer[0] != 'G' && kmer[0] != 'T') continue;
        }
        // Pack the k-mer into 10 bytes.
        unsigned char kbytes[10];
        for (int i = 0; i < 10; i++) {
            unsigned int b = 0;
            for (int k = 0; k < 4; k++) {
                b = (b << 2) | base_to_bits[(unsigned char)kmer[i*4+k]];
            }
            kbytes[i] = (unsigned char)b;
        }
        // Tokenize the rest: skip kmer, skip mask, skip lcp, read sign, contig, '|', position.
        char *q = kmer + 40;
        // mask
        while (*q == ' ') q++;
        while (*q && *q != ' ' && *q != '\n') q++;
        // lcp
        while (*q == ' ') q++;
        while (*q && *q != ' ' && *q != '\n') q++;
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
        int strand = (sign == '+') ? 0 : 1;
        // Build output: "contig position strand hex"
        int n = snprintf(out, sizeof(out), "%s %s %d ", contig_start, pos_start, strand);
        for (int i = 0; i < 10; i++) {
            out[n++] = hex_chars[kbytes[i] >> 4];
            out[n++] = hex_chars[kbytes[i] & 0xf];
        }
        out[n++] = '\n';
        fwrite(out, 1, n, stdout);
    }
    return 0;
}
