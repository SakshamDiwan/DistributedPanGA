// pga-extract: Phase 1 single-process syncmer extractor.
// Usage:  pga-extract <gdb-stem>  [output-tuples-file]
// Output format (one line per tuple, sorted):
//   <contig_id> <position> <strand> <20-hex-character k-mer>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "extract.h"
#include "constants.h"
#include "GDB.h"

static const char hex_chars[] = "0123456789abcdef";

static void dump_tuples(FILE *out, const TupleBuffer *buf)
{
    char line[128];
    for (int64_t i = 0; i < buf->count; i++)
    {
        const Tuple *t = &buf->data[i];
        int n = snprintf(line, sizeof(line), "%d %lld %d ",
                         t->contig, (long long) t->position, t->strand);
        for (int b = 0; b < KBYTES; b++)
        {
            line[n++] = hex_chars[t->kmer[b] >> 4];
            line[n++] = hex_chars[t->kmer[b] & 0xf];
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

    TupleBuffer buf;
    tuple_buffer_init(&buf);
    extract_to_tuples(gdb, /*num_workers=*/1, /*worker_id=*/0, &buf);

    fprintf(stderr, "Extracted %lld tuples\n", (long long) buf.count);
    // Output is intentionally NOT sorted here. Pipe through `sort` externally
    // to match gold.tuples' lexicographic text-sort order.

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
    dump_tuples(out, &buf);

    if (out != stdout) fclose(out);
    tuple_buffer_free(&buf);
    Close_GDB(gdb);
    return 0;
}
