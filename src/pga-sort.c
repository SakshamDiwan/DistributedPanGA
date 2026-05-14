// pga-sort: Phase 2 single-process syncmer extractor + msd_sort.
// Usage:  pga-sort <gdb-stem>  [output-tuples-file]
//
// Pipeline:
//   1. Read GDB.
//   2. Extract records into FastGA-compatible byte layout (RecordBuffer).
//   3. msd_sort the buffer (sorts by k-mer key, fills LCP byte).
//   4. Decode each record back to text tuple format and dump.
//
// Output format (one line per record, sort order is k-mer-byte order):
//   <contig_id> <position> <strand> <20-hex-character k-mer>
//
// To verify byte-equivalence with FastGA's index, sort the output and
// compare against gold.tuples (exactly as for Phase 1).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "extract.h"
#include "record.h"
#include "sort.h"
#include "constants.h"
#include "GDB.h"

static const char hex_chars[] = "0123456789abcdef";

static int64_t read_le_field(const uint8_t *bytes, int n)
{
    int64_t v = 0;
    for (int k = 0; k < n; k++)
        v |= ((int64_t) bytes[k]) << (8 * k);
    return v;
}

static void dump_records_as_tuples(FILE *out, const RecordBuffer *buf)
{
    int rsize       = buf->sizing.record_size;
    int post_bytes  = buf->sizing.post_bytes;
    int cont_bytes  = buf->sizing.cont_bytes;
    int pos_off     = POS_OFFSET;
    int cont_off    = POS_OFFSET + post_bytes;

    // Strand bit lives in the high bit of the topmost contig byte.
    int64_t strand_mask = ((int64_t) 1) << (8 * cont_bytes - 1);
    int64_t contig_mask = strand_mask - 1;

    char line[128];
    for (int64_t i = 0; i < buf->count; i++)
    {
        const uint8_t *r = buf->data + i * rsize;

        int64_t pos     = read_le_field(r + pos_off, post_bytes);
        int64_t cont_v  = read_le_field(r + cont_off, cont_bytes);
        int     strand  = (cont_v & strand_mask) ? 1 : 0;
        int     contig  = (int) (cont_v & contig_mask);

        int n = snprintf(line, sizeof(line), "%d %lld %d ",
                         contig, (long long) pos, strand);
        for (int b = 0; b < KBYTES; b++)
        {
            uint8_t kb = r[KMER_OFFSET + b];
            line[n++] = hex_chars[kb >> 4];
            line[n++] = hex_chars[kb & 0xf];
        }
        line[n++] = '\n';
        fwrite(line, 1, n, out);
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2 || argc > 3)
    {
        fprintf(stderr, "Usage: %s <gdb-stem> [output-tuples-file]\n", argv[0]);
        return 1;
    }

    GDB _gdb, *gdb = &_gdb;
    if (Read_GDB(gdb, argv[1]) < 0)
    {
        fprintf(stderr, "Read_GDB failed for '%s'\n", argv[1]);
        return 1;
    }

    fprintf(stderr, "Loaded GDB: %d contigs, %lld total bases\n",
            gdb->ncontig, (long long) gdb->seqtot);

    RecordSizing sizing = compute_record_sizing(gdb);
    fprintf(stderr, "Record layout: post_bytes=%d cont_bytes=%d record_size=%d\n",
            sizing.post_bytes, sizing.cont_bytes, sizing.record_size);

    RecordBuffer buf;
    record_buffer_init(&buf, sizing);

    fprintf(stderr, "Extracting records...\n");
    extract_to_records(gdb, /*num_workers=*/1, /*worker_id=*/0, &buf);
    fprintf(stderr, "Extracted %lld records\n", (long long) buf.count);

    fprintf(stderr, "Sorting via msd_sort...\n");
    sort_records(&buf);

    // Verify msd_sort actually produced k-mer-byte-lex order on our record layout.
    // Walks the buffer in place; tie on full k-mer is allowed (memcmp == 0).
    {
        int rsize = buf.sizing.record_size;
        int64_t bad = 0;
        int64_t first_bad = -1;
        for (int64_t i = 1; i < buf.count; i++)
        {
            const uint8_t *a = buf.data + (i - 1) * rsize + KMER_OFFSET;
            const uint8_t *b = buf.data +  i      * rsize + KMER_OFFSET;
            if (memcmp(a, b, KBYTES) > 0)
            {
                if (first_bad < 0) first_bad = i;
                bad++;
            }
        }
        if (bad == 0)
            fprintf(stderr, "Sort order check: PASS (%lld records non-decreasing)\n",
                    (long long) buf.count);
        else
        {
            fprintf(stderr, "Sort order check: FAIL (%lld violations, first at index %lld)\n",
                    (long long) bad, (long long) first_bad);
            return 2;
        }
    }


    FILE *out = stdout;
    if (argc == 3)
    {
        out = fopen(argv[2], "w");
        if (out == NULL)
        {
            fprintf(stderr, "Cannot open '%s' for write\n", argv[2]);
            return 1;
        }
    }

    fprintf(stderr, "Writing tuples...\n");
    dump_records_as_tuples(out, &buf);

    if (out != stdout) fclose(out);
    record_buffer_free(&buf);
    Close_GDB(gdb);
    return 0;
}
