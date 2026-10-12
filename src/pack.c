#include "pack.h"
#include "constants.h"
#include "tables.h"
#include "extract.h"

int pack_forward_kmer(GDB *gdb, int contig_id, int64_t j, uint8_t *out)
{
    int64_t begb = j;
    int64_t endb = j + KMER;          // exclusive
    int     clen = gdb->contigs[contig_id].clen;
    // Get_Contig_Piece writes a sentinel at buffer[-1], so we pass buf+1.
    char    raw[65];
    char   *buf = raw + 1;
    uint8_t *bases;

    if (begb < 0 || endb > clen)
        return -1;

    g_extract_stats.pack_reads += 1;
    g_extract_stats.pack_bases += endb - begb;
    bases = (uint8_t *) Get_Contig_Piece(gdb, contig_id, (int) begb, (int) endb, NUMERIC, buf);
    if (bases == NULL)
        return -1;

    pack_forward_bases(bases, out);
    return 0;
}

void pack_forward_bases(const uint8_t *bases, uint8_t *out)
{
    int i;

    for (i = 0; i < KBYTES; i++)
    {
        int p = i * 4;
        out[i] = (uint8_t) ((bases[p]   << 6) | (bases[p+1] << 4)
                          | (bases[p+2] << 2) |  bases[p+3]);
    }
}

int pack_rc_kmer(GDB *gdb, int contig_id, int64_t j, uint8_t *out)
{
    int64_t begb = j - 28;
    int64_t endb = j + 12;            // exclusive
    int     clen = gdb->contigs[contig_id].clen;
    // Get_Contig_Piece writes a sentinel at buffer[-1], so we pass buf+1.
    char    raw[65];
    char   *buf = raw + 1;
    uint8_t *bases;

    if (begb < 0 || endb > clen)
        return -1;

    g_extract_stats.pack_reads += 1;
    g_extract_stats.pack_bases += endb - begb;
    bases = (uint8_t *) Get_Contig_Piece(gdb, contig_id, (int) begb, (int) endb, NUMERIC, buf);
    if (bases == NULL)
        return -1;

    pack_rc_bases(bases, out);
    return 0;
}

void pack_rc_bases(const uint8_t *bases, uint8_t *out)
{
    int b;

    // bases[k] = forward NUMERIC base at contig position (j-28 + k), k in [0, 40).
    // For byte b in [0, KBYTES), the four source positions are [j+8-4b, j+11-4b],
    // which become local indices [36-4b, 39-4b].
    for (b = 0; b < KBYTES; b++)
    {
        int p3 = 36 - 4*b;
        uint8_t fwd4 = (uint8_t) ((bases[p3]   << 6) | (bases[p3+1] << 4)
                                | (bases[p3+2] << 2) |  bases[p3+3]);
        out[b] = (uint8_t) Comp[fwd4];
    }
}
