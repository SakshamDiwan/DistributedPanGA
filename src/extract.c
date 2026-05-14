#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "extract.h"
#include "constants.h"
#include "tables.h"
#include "pack.h"
#include "work_split.h"
#include "GDB.h"

#define SCAN_MAX 10000000   // must be divisible by 4

void tuple_buffer_init(TupleBuffer *b)
{
    b->capacity = 1 << 20;
    b->count    = 0;
    b->data     = (Tuple *) malloc(b->capacity * sizeof(Tuple));
    if (b->data == NULL)
    {
        fprintf(stderr, "tuple_buffer_init: out of memory\n");
        exit(1);
    }
}

void tuple_buffer_free(TupleBuffer *b)
{
    free(b->data);
    b->data = NULL;
    b->count = 0;
    b->capacity = 0;
}

static void tuple_buffer_grow(TupleBuffer *b)
{
    b->capacity *= 2;
    b->data = (Tuple *) realloc(b->data, b->capacity * sizeof(Tuple));
    if (b->data == NULL)
    {
        fprintf(stderr, "tuple_buffer_grow: out of memory\n");
        exit(1);
    }
}

static inline void emit_tuple(TupleBuffer *buf, int contig, int64_t position, int strand,
                              const uint8_t *kmer)
{
    if (buf->count >= buf->capacity)
        tuple_buffer_grow(buf);
    Tuple *t = &buf->data[buf->count++];
    t->contig   = contig;
    t->position = position;
    t->strand   = strand;
    memcpy(t->kmer, kmer, KBYTES);
}

// Generic emit callback: invoked once per syncmer emission with the validated
// packed k-mer and the stored position.  strand=0 for forward, 1 for RC.
typedef void (*EmitFn)(void *user, int contig, int64_t position, int strand,
                       const uint8_t *kmer);

// Per-contig scan. Mirrors scan_thread in notes/GIXmake_dump.c (line 486+).
// Variable names and control flow are kept verbatim where possible.
static void scan_contig(GDB *gdb, int r, char *SEQ_buf, EmitFn emit, void *user)
{
    int len, beg, end;
    int i, n;
    uint8_t *seq, *sn;
    int ne[8];
    int nq[8], cq[8];
    int mzr[8];
    int min4, pos4;
    int lnK, KMT;

    if (gdb->contigs[r].boff < 0)
        return;

    len = gdb->contigs[r].clen;
    if (len > SCAN_MAX)
        end = SCAN_MAX;
    else
        end = len;

    seq = (uint8_t *) Get_Contig_Piece(gdb, r, 0, end, NUMERIC, SEQ_buf);

    lnK = len - KMER;
    KMT = KMER - TMER;

    sn = seq + 3;
    n = (seq[0] << 4) | (seq[1] << 2) | seq[2];
    for (i = 0; i < 4; i++)
    {
        int c;
        ne[i+4] = n = ((n << 2) | sn[i]) & 0xff;
        c = Comp[n];
        nq[i] = TMap[n];
        cq[i] = TMap[c];
    }
    sn += 4;
    min4 = 0x10000;
    pos4 = 0;
    for (i = 0; i < SOFF; i++)
    {
        int nh, ch, c;
        int mn, mc, mz;

        ne[i] = n = ((n << 2) | sn[i]) & 0xff;
        c  = Comp[n];
        nh = TMap[n];
        ch = TMap[c];
        mn = (nq[i] << 8) | nh;
        mc = cq[i] | (ch << 8);
        nq[i] = nh;
        cq[i] = ch;
        if (mn < mc)
            mzr[i] = mz = mn;
        else
            mzr[i] = mz = mc;
        if (mz < min4)
        {
            min4 = mz;
            pos4 = i;
        }
    }

    beg = SOFF;
    while (1)
    {
        end -= SMER;
        for (i = beg; i <= end; i++)
        {
            int iq, jq;
            int nh, ch;
            int mn, mc, mz;
            int j, w, c, p;

            iq = i & 0x7;
            w = ne[iq];
            ne[iq] = n = ((n << 2) | sn[i]) & 0xff;
            c  = Comp[n];
            nh = TMap[n];
            ch = TMap[c];

            jq = i & 0x3;
            mn = (nq[jq] << 8) | nh;
            mc = cq[jq] | (ch << 8);
            nq[jq] = nh;
            cq[jq] = ch;
            if (mn < mc)
                mz = mzr[jq] = mn;
            else
                mz = mzr[jq] = mc;

            if (mz < min4)                          // right-end of 12-syncmer
            {
                min4 = mz;
                pos4 = i;
            }
            else if (pos4 == i - SOFF)              // left-end of 12-syncmer
            {
                min4 = mzr[(++pos4) & 0x3];
                for (j = pos4 + 1; j <= i; j++)
                    if (mzr[j & 0x3] < min4)
                    {
                        min4 = mzr[j & 0x3];
                        pos4 = j;
                    }
            }
            else if (mz > min4)
                continue;
            // mz == min4 falls through (the "Hit RE" tie branch)

            j = i - SOFF;
            p = ne[(i+4) & 0x7];
            (void) w;  // w isn't needed here; the bucket is recomputed in Phase 2 from packed bytes
            (void) c;  // same reason

            if (j <= lnK)
            {
                uint8_t kmer_buf[KBYTES];
                if (pack_forward_kmer(gdb, r, j, kmer_buf) == 0)
                    emit(user, r, j, 0, kmer_buf);
            }
            if (j >= KMT)
            {
                uint8_t kmer_buf[KBYTES];
                if (pack_rc_kmer(gdb, r, j, kmer_buf) == 0)
                    emit(user, r, j + TMER, 1, kmer_buf);
            }
        }
        end += SMER;

        if (end == len)
            break;
        beg = end;
        end += SCAN_MAX;
        if (end > len)
            end = len;
        beg -= (SMER - 1);
        sn = (uint8_t *) Get_Contig_Piece(gdb, r,
                                          beg + (SMER - 1), end,
                                          NUMERIC, SEQ_buf) - beg;
    }
}

// --- TupleBuffer emit path (Phase 1) ---------------------------------------

static void tuple_emit_cb(void *user, int contig, int64_t position, int strand,
                          const uint8_t *kmer)
{
    emit_tuple((TupleBuffer *) user, contig, position, strand, kmer);
}

// --- RecordBuffer emit path (Phase 2) --------------------------------------

void record_buffer_init(RecordBuffer *b, RecordSizing sizing)
{
    b->sizing   = sizing;
    b->capacity = 1 << 20;
    b->count    = 0;
    // +1 for msd_sort's sentinel write at array[asize].
    b->data = (uint8_t *) malloc(b->capacity * sizing.record_size + 1);
    if (b->data == NULL)
    {
        fprintf(stderr, "record_buffer_init: out of memory\n");
        exit(1);
    }
}

void record_buffer_free(RecordBuffer *b)
{
    free(b->data);
    b->data = NULL;
    b->count = 0;
    b->capacity = 0;
}

static void record_buffer_grow(RecordBuffer *b)
{
    b->capacity *= 2;
    b->data = (uint8_t *) realloc(b->data, b->capacity * b->sizing.record_size + 1);
    if (b->data == NULL)
    {
        fprintf(stderr, "record_buffer_grow: out of memory\n");
        exit(1);
    }
}

static void record_emit_cb(void *user, int contig, int64_t position, int strand,
                           const uint8_t *kmer)
{
    RecordBuffer *b = (RecordBuffer *) user;
    if (b->count >= b->capacity)
        record_buffer_grow(b);

    int post_bytes  = b->sizing.post_bytes;
    int cont_bytes  = b->sizing.cont_bytes;
    int record_size = b->sizing.record_size;
    uint8_t *r = b->data + b->count * record_size;

    r[LCP_OFFSET] = 0;                 // LCP placeholder
    memcpy(r + KMER_OFFSET, kmer, KBYTES);
    r[MASK_OFFSET] = 0;                // mask: not supported in v1

    int64_t pos_val = position;
    for (int k = 0; k < post_bytes; k++)
        r[POS_OFFSET + k] = (uint8_t) ((pos_val >> (8 * k)) & 0xff);

    // Contig field: contig_id with strand bit in the high bit of the last cont_byte.
    int64_t cont_val = (int64_t) contig;
    if (strand)
        cont_val |= ((int64_t) 1) << (8 * cont_bytes - 1);
    int cont_off = POS_OFFSET + post_bytes;
    for (int k = 0; k < cont_bytes; k++)
        r[cont_off + k] = (uint8_t) ((cont_val >> (8 * k)) & 0xff);

    b->count++;
}

// --- Driver functions ------------------------------------------------------

static void run_scan(GDB *gdb, int num_workers, int worker_id, EmitFn emit, void *user)
{
    int *split = (int *) malloc((num_workers + 1) * sizeof(int));
    int64_t *post = (int64_t *) malloc((num_workers + 1) * sizeof(int64_t));
    compute_work_split(gdb, num_workers, split, post);

    init_comp_table();

    char *SEQ_buf = (char *) malloc(SCAN_MAX + 8);
    if (SEQ_buf == NULL)
    {
        fprintf(stderr, "run_scan: out of memory for SEQ buffer\n");
        exit(1);
    }
    char *SEQ = SEQ_buf + 1;        // Get_Contig_Piece writes sentinel at buffer[-1]

    int rbeg = split[worker_id];
    int rend = split[worker_id + 1];
    for (int r = rbeg; r < rend; r++)
        scan_contig(gdb, r, SEQ, emit, user);

    free(SEQ_buf);
    free(split);
    free(post);
}

void extract_to_tuples(GDB *gdb, int num_workers, int worker_id, TupleBuffer *buf)
{
    run_scan(gdb, num_workers, worker_id, tuple_emit_cb, buf);
}

void extract_to_records(GDB *gdb, int num_workers, int worker_id, RecordBuffer *buf)
{
    run_scan(gdb, num_workers, worker_id, record_emit_cb, buf);
}
