#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ingest.h"
#include "gene_core.h"

// Syncmer parameters (same as FastGA)
#define TMER 12
#define SMER  8
#define SOFF  4   // TMER - SMER

// DNA complement table (byte level)
static int Comp[256];

// Syncmer hash table (same as FastGA's TMap)
static int TMap[256] =
 { 0xff, 0xd4, 0xf5, 0xfd, 0xe4, 0xad, 0x21, 0xa5, 0xed, 0x64, 0xbf, 0xa9, 0xf3, 0x70, 0xd6, 0xf0,
   0xca, 0x89, 0xcb, 0xc9, 0x82, 0x9d, 0x13, 0x79, 0x0a, 0x0f, 0x25, 0x19, 0x3e, 0x47, 0xa3, 0xa8,
   0xf9, 0x5e, 0xe8, 0xa1, 0xb0, 0x71, 0x1d, 0x8c, 0xde, 0x69, 0xe7, 0x7c, 0x56, 0x3f, 0x90, 0xa4,
   0xeb, 0x45, 0x59, 0xf1, 0x97, 0x4c, 0x08, 0xa0, 0xb8, 0x4a, 0x86, 0xc8, 0xcd, 0x98, 0x7d, 0xfc,
   0xef, 0x4d, 0x83, 0x7e, 0xdc, 0x66, 0x2b, 0x8e, 0xe0, 0xa7, 0xd0, 0xa2, 0x88, 0x5f, 0x7f, 0xd9,
   0x9b, 0x78, 0xd1, 0x8b, 0xc3, 0x8f, 0x2d, 0xe6, 0x18, 0x27, 0x2c, 0x24, 0x94, 0xb7, 0xce, 0xbd,
   0x0d, 0x04, 0x1c, 0x09, 0x16, 0x23, 0x00, 0x1e, 0x1a, 0x29, 0x2e, 0x15, 0x01, 0x10, 0x2a, 0x20,
   0xbe, 0x31, 0x43, 0x58, 0xc2, 0xaa, 0x1f, 0xe5, 0xc5, 0x9e, 0xcf, 0xc6, 0x68, 0xb2, 0x80, 0xf4,
   0xf8, 0x53, 0xb6, 0x93, 0x76, 0x37, 0x11, 0x40, 0xda, 0x51, 0xba, 0x46, 0x42, 0x30, 0x60, 0x6d,
   0x5c, 0x39, 0x9f, 0x48, 0x6c, 0x62, 0x28, 0x67, 0x06, 0x12, 0x26, 0x0e, 0x33, 0x50, 0xa6, 0x63,
   0xdd, 0x3b, 0xab, 0x4b, 0x72, 0x5b, 0x22, 0x6f, 0xb4, 0x61, 0x92, 0x99, 0x36, 0x38, 0x65, 0xac,
   0x4f, 0x2f, 0x32, 0x44, 0x54, 0x3c, 0x03, 0x5d, 0x73, 0x3a, 0x77, 0x84, 0x8d, 0x4e, 0x49, 0xd2,
   0xfb, 0x91, 0x6a, 0xcc, 0x8a, 0x35, 0x02, 0x55, 0x7a, 0x34, 0x96, 0x3d, 0xd3, 0x41, 0x85, 0xf2,
   0xb1, 0x75, 0xc4, 0xb5, 0xbb, 0xb3, 0x1b, 0xd5, 0x07, 0x05, 0x17, 0x0b, 0x7b, 0xd7, 0xdf, 0xea,
   0xe3, 0x57, 0xc0, 0x95, 0x9c, 0x6e, 0x14, 0xae, 0xb9, 0x6b, 0xc1, 0x81, 0x87, 0x74, 0xd8, 0xe2,
   0xec, 0x52, 0xbc, 0xe9, 0xe1, 0xdb, 0x0c, 0xf7, 0xaf, 0x5a, 0x9a, 0xc7, 0xfa, 0xf6, 0xee, 0xfe
 };

// Initialize complement table (same logic as GIXmake)
static void init_complement_table(void)
{
    int i, l0, l1, l2, l3;
    i = 0;
    for (l0 = 3; l0 >= 0; l0 -= 1)
     for (l1 = 12; l1 >= 0; l1 -= 4)
      for (l2 = 48; l2 >= 0; l2 -= 16)
       for (l3 = 192; l3 >= 0; l3 -= 64)
         Comp[i++] = (l3 | l2 | l1 | l0);
}

// Compute how many bytes needed to represent a value
static int bytes_needed(int64_t max_val)
{
    int bytes = 0;
    int64_t cum = 1;
    while (cum < max_val)
    {   cum *= 256;
        bytes += 1;
    }
    return bytes;
}

void kmer_buffer_init(KmerBuffer *buf, int record_size, int64_t initial_capacity)
{
    buf->record_size = record_size;
    buf->capacity = initial_capacity;
    buf->count = 0;
    buf->records = (uint8_t *) malloc(initial_capacity * record_size);
    if (buf->records == NULL)
    {   fprintf(stderr, "Failed to allocate k-mer buffer\n");
        exit(1);
    }
}

void kmer_buffer_free(KmerBuffer *buf)
{
    free(buf->records);
    buf->records = NULL;
    buf->count = 0;
    buf->capacity = 0;
}

// Double the buffer if full
static void kmer_buffer_grow(KmerBuffer *buf)
{
    buf->capacity *= 2;
    buf->records = (uint8_t *) realloc(buf->records, buf->capacity * buf->record_size);
    if (buf->records == NULL)
    {   fprintf(stderr, "Failed to grow k-mer buffer\n");
        exit(1);
    }
}

// Compute contig split across MPI ranks (same logic as GIXmake's DBsplit)
static void compute_split(GDB *gdb, int num_ranks, int *split, int64_t *post)
{
    int64_t total_bases = gdb->seqtot;
    int64_t bases_so_far = 0;
    int64_t next_threshold;
    int next_rank = 1;
    int r;

    split[0] = 0;
    post[0]  = 0;
    next_threshold = total_bases / num_ranks;

    for (r = 0; r < gdb->ncontig; r++)
    {   bases_so_far += gdb->contigs[r].clen;
        while (bases_so_far >= next_threshold && next_rank < num_ranks)
        {   split[next_rank] = r + 1;
            post[next_rank]  = bases_so_far;
            next_rank += 1;
            next_threshold = (total_bases * next_rank) / num_ranks;
        }
    }
    split[num_ranks] = gdb->ncontig;
    post[num_ranks]  = total_bases;
}

void extract_kmers(GDB *gdb, int rank, int num_ranks, int kmer, KmerBuffer *buf, int64_t *counts)
{
    int kbytes = kmer / 4;
    int post_bytes = bytes_needed(gdb->maxctg); 
    int cont_bytes = bytes_needed(2 * gdb->ncontig);
    int record_size = 1 + kbytes + 1 + post_bytes + cont_bytes;  // LCP + kmer + mask + pos + contig

    // Compute work split
    int     *split = (int *) malloc((num_ranks + 1) * sizeof(int));
    int64_t *posts = (int64_t *) malloc((num_ranks + 1) * sizeof(int64_t));
    compute_split(gdb, num_ranks, split, posts);

    printf("Rank %d: assigned contigs [%d, %d)\n", rank, split[rank], split[rank+1]);
    fflush(stdout);

    // Init complement table
    init_complement_table();

    // Init output buffer
    kmer_buffer_init(buf, record_size, 1000000);

    memset(counts, 0, 1024 * sizeof(int64_t));

    // Scan assigned contigs
    int seq_max = 10000000;
    char *seq_buf = (char *) malloc(seq_max + 8);

    for (int r = split[rank]; r < split[rank+1]; r++)
    {
        int len = gdb->contigs[r].clen;
        if (gdb->contigs[r].boff < 0)
            continue;   // skip fake contigs from short_GDB_fix

        int end = (len > seq_max) ? seq_max : len;
        uint8_t *seq = (uint8_t *) Get_Contig_Piece(gdb, r, 0, end, NUMERIC, seq_buf + 1);

        // Pack initial bytes
        uint8_t *bases = seq + 3;
        int byte_window = (seq[0] << 4) | (seq[1] << 2) | seq[2];

        int window_history[8], fwd_hash[8], rev_hash[8], syncmer_score[8];
        int best_score, best_pos;

        for (int i = 0; i < 4; i++)
        {   window_history[i+4] = byte_window = ((byte_window << 2) | bases[i]) & 0xff;
            int rev_byte = Comp[byte_window];
            fwd_hash[i] = TMap[byte_window];
            rev_hash[i] = TMap[rev_byte];
        }
        bases += 4;
        best_score = 0x10000;
        best_pos = 0;

        for (int i = 0; i < SOFF; i++)
        {   window_history[i] = byte_window = ((byte_window << 2) | bases[i]) & 0xff;
            int rev_byte     = Comp[byte_window];
            int fwd_score    = TMap[byte_window];
            int rev_score    = TMap[rev_byte];
            int fwd_syncmer  = (fwd_hash[i] << 8) | fwd_score;
            int rev_syncmer  = rev_hash[i] | (rev_score << 8);
            fwd_hash[i] = fwd_score;
            rev_hash[i] = rev_score;
            int min_syncmer;
            if (fwd_syncmer < rev_syncmer)
                syncmer_score[i] = min_syncmer = fwd_syncmer;
            else
                syncmer_score[i] = min_syncmer = rev_syncmer;
            if (min_syncmer < best_score)
            {   best_score = min_syncmer;
                best_pos = i;
            }
        }

        int last_kmer_pos = len - kmer;
        int KMT = kmer - TMER;
        int beg = SOFF;

        while (1)
        {   end -= SMER;
            for (int i = beg; i <= end; i++)
            {
                int iq = i & 0x7;
                int w  = window_history[iq];
                window_history[iq] = byte_window = ((byte_window << 2) | bases[i]) & 0xff;
                int rev_byte  = Comp[byte_window];
                int fwd_score = TMap[byte_window];
                int rev_score = TMap[rev_byte];

                int jq = i & 0x3;
                int fwd_syncmer = (fwd_hash[jq] << 8) | fwd_score;
                int rev_syncmer = rev_hash[jq] | (rev_score << 8);
                fwd_hash[jq] = fwd_score;
                rev_hash[jq] = rev_score;
                int min_syncmer;
                if (fwd_syncmer < rev_syncmer)
                    min_syncmer = syncmer_score[jq] = fwd_syncmer;
                else
                    min_syncmer = syncmer_score[jq] = rev_syncmer;

                if (min_syncmer < best_score)
                {   best_score = min_syncmer;
                    best_pos = i;
                }
                else if (best_pos == i - SOFF)
                {   best_score = syncmer_score[(++best_pos) & 0x3];
                    for (int j = best_pos + 1; j <= i; j++)
                        if (syncmer_score[j & 0x3] < best_score)
                        {   best_score = syncmer_score[j & 0x3];
                            best_pos = j;
                        }
                }
                else if (min_syncmer > best_score)
                    continue;

                // Found a syncmer at position j = i - SOFF
                int j = i - SOFF;
                int p = window_history[(i+4) & 0x7];

                // Pack forward strand k-mer
                if (j <= last_kmer_pos)
                {
                    if (buf->count >= buf->capacity)
                        kmer_buffer_grow(buf);

                    uint8_t *x = buf->records + buf->count * buf->record_size;

                    *x++ = 0;  // LCP placeholder at byte 0 (msd_sort writes LCP here)
                    for (int k = 0; k < kbytes; k++)
                    {   int base = j + k*4;
                        *x++ = (seq[base]<<6) | (seq[base+1]<<4) | (seq[base+2]<<2) | seq[base+3];
                    }
                    *x++ = 0;  // mask byte: always 0 — we do not read .1ano repeat-mask annotations.
                               // FastGA sets this to the number of leading bases of the k-mer that
                               // fall inside a masked region; the merge phase skips seeds where this
                               // value >= the current LCP depth. For yeast (few repeats) this is fine.
                               // TODO: add mask support for large repeat-rich genomes (e.g. human).
                               // Without masking, we will produce more seeds than FastGA in those cases.

                    // Position within contig
                    int64_t pos_val = j;
                    for (int k = 0; k < post_bytes; k++)
                        *x++ = ((uint8_t *)&pos_val)[k];

                    // Contig ID (no strand bit for forward)
                    int64_t cont_val = r;
                    for (int k = 0; k < cont_bytes; k++)
                        *x++ = ((uint8_t *)&cont_val)[k];

                    counts[w << 2 | p >> 6]++;
                    buf->count++;
                }

                // Pack reverse complement k-mer
                if (j >= KMT)
                {
                    if (buf->count >= buf->capacity)
                        kmer_buffer_grow(buf);

                    uint8_t *x = buf->records + buf->count * buf->record_size;

                    *x++ = 0;  // LCP placeholder at byte 0 (msd_sort writes LCP here)
                    int rc_start = j + TMER;
                    for (int k = 0; k < kbytes; k++)
                    {   int base = rc_start - k*4;
                        *x++ = Comp[(seq[base-3]<<6) | (seq[base-2]<<4) | (seq[base-1]<<2) | seq[base]];
                    }
                    *x++ = 0;  // mask byte: always 0, see note in forward strand branch above.

                    int64_t pos_val = j;
                    for (int k = 0; k < post_bytes; k++)
                        *x++ = ((uint8_t *)&pos_val)[k];

                    int64_t cont_val = (int64_t)r | ((int64_t)1 << (cont_bytes*8 - 1));
                    for (int k = 0; k < cont_bytes; k++)
                        *x++ = ((uint8_t *)&cont_val)[k];

                    counts[rev_byte << 2 | (3 - (p & 0x3))]++;
                    buf->count++;
                }
            }
            end += SMER;

            if (end == len)
                break;
            beg = end;
            end += seq_max;
            if (end > len)
                end = len;
            beg -= (SMER - 1);
            bases = (uint8_t *) Get_Contig_Piece(gdb, r, beg + (SMER-1), end, NUMERIC, seq_buf + 1) - beg;
        }
    }

    printf("Rank %d: extracted %lld k-mers from %d contigs\n",
           rank, buf->count, split[rank+1] - split[rank]);
    fflush(stdout);

    free(seq_buf);
    free(split);
    free(posts);
}