#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kmer_adapter.h"

int build_kmer_adapter(const RecordBuffer *buf,
                       KmerStreamAdapter *T,
                       PostListAdapter *P)
{
    int64 nels       = buf->count;
    int   rec_size   = buf->sizing.record_size;
    int   post_bytes = buf->sizing.post_bytes;
    int   cont_bytes = buf->sizing.cont_bytes;

    T->kmer   = 40;
    T->minval = 1;
    T->ibyte  = 3;
    T->kbyte  = 10;
    T->hbyte  = T->kbyte - T->ibyte;                              // 7
    T->pbyte  = T->hbyte + 1 + 1 + post_bytes + cont_bytes;       // 13
    T->tbyte  = T->ibyte + T->pbyte;                              // 16, used by Current_Entry
    T->nels   = nels;
    T->ixlen  = 1 << (8 * T->ibyte);
    T->shift  = 0;
    T->inver  = NULL;
    T->copn   = -1;
    T->part   = 0;
    T->nthr   = 1;
    T->nlen   = 0;
    T->name   = NULL;
    T->neps   = NULL;

    T->table = (uint8 *) malloc(nels * T->pbyte);
    if (T->table == NULL) {
        fprintf(stderr, "build_kmer_adapter: out of memory for table\n");
        return -1;
    }

    T->index = (int64 *) calloc(T->ixlen, sizeof(int64));
    if (T->index == NULL) {
        fprintf(stderr, "build_kmer_adapter: out of memory for index\n");
        free(T->table);
        return -1;
    }

    uint8 *dst = T->table;
    int64 max_run = 0, current_run = 0;
    int   prev_prefix = -1;

    for (int64 i = 0; i < nels; i++) {
        const uint8 *src = buf->data + i * rec_size;

        int prefix = (src[KMER_OFFSET] << 16) |
                     (src[KMER_OFFSET + 1] << 8) |
                     src[KMER_OFFSET + 2];

        if (prefix == prev_prefix) {
            current_run++;
        } else {
            if (current_run > max_run) max_run = current_run;
            current_run = 1;
            prev_prefix = prefix;
        }

        memcpy(dst, src + KMER_OFFSET + T->ibyte, T->hbyte);
        dst += T->hbyte;
        *dst++ = 0;
        *dst++ = src[LCP_OFFSET];
        memcpy(dst, src + POS_OFFSET, post_bytes + cont_bytes);
        dst += post_bytes + cont_bytes;
    }
    if (current_run > max_run) max_run = current_run;

    // "≤ p" cumulative: index[p] = number of records with prefix ≤ p.
    // Matches libfastk.c semantics (see First_Kmer_Entry, Next_Kmer_Entry).
    memset(T->index, 0, T->ixlen * sizeof(int64));
    for (int64 i = 0; i < nels; i++) {
        const uint8 *src = buf->data + i * rec_size;
        int prefix = (src[KMER_OFFSET] << 16) |
                    (src[KMER_OFFSET + 1] << 8) |
                    src[KMER_OFFSET + 2];
        T->index[prefix]++;
    }
    {
        int64 cum = 0;
        for (int64 p = 0; p < T->ixlen; p++) {
            cum += T->index[p];
            T->index[p] = cum;
        }
    }

    T->csuf  = T->table;
    T->cidx  = 0;
    T->cpre  = 0;
    while (T->cpre < T->ixlen && T->index[T->cpre] <= 0)
        T->cpre++;
    T->clone = 0;
    T->ctop = T->table + nels * T->pbyte;
    P->pbyte = post_bytes + cont_bytes;
    P->cbyte = cont_bytes;
    P->nels  = nels;
    P->maxp  = max_run + 1;
    P->freq  = 0;
    P->nctg  = 0;
    P->perm  = NULL;

    return 0;
}

void free_kmer_adapter(KmerStreamAdapter *T, PostListAdapter *P)
{
    (void) P;
    if (T->table) { free(T->table); T->table = NULL; }
    if (T->index) { free(T->index); T->index = NULL; }
    // T->inver is NULL — kept for layout compat with Kmer_Stream
}

void adapter_First_Kmer_Entry(KmerStreamAdapter *T)
{
    T->csuf = T->table;
    T->cidx = 0;
    T->cpre = 0;
    while (T->cpre < T->ixlen && T->index[T->cpre] <= 0)
        T->cpre++;
}

void adapter_Next_Kmer_Entry(KmerStreamAdapter *T)
{
    T->csuf += T->pbyte;
    T->cidx++;
    if (T->csuf >= T->ctop) {
        T->csuf = NULL;
        T->cpre = T->ixlen;
        return;
    }
    while (T->cpre < T->ixlen && T->index[T->cpre] <= T->cidx)
        T->cpre++;
}

// O(ixlen) linear scan — acceptable only because self_adaptamer_merge
// only calls this NTHREADS-1 times. If NTHREADS grows, build T->inver
// in build_kmer_adapter and use S->inver[i>>S->shift] like libfastk.c:1289.
void adapter_GoTo_Kmer_Index(KmerStreamAdapter *T, int64 i)
{
    if (i >= T->nels) {
        T->csuf = NULL;
        T->cpre = T->ixlen;
        T->cidx = T->nels;
        return;
    }
    T->cidx = i;
    T->csuf = T->table + i * T->pbyte;
    T->cpre = 0;
    while (T->cpre < T->ixlen && T->index[T->cpre] <= i)
        T->cpre++;
}

KmerStreamAdapter *adapter_Clone_Kmer_Stream(KmerStreamAdapter *T)
{
    KmerStreamAdapter *C = malloc(sizeof(KmerStreamAdapter));
    *C = *T;
    C->clone = 1;
    return C;
}

void adapter_Free_Kmer_Stream(KmerStreamAdapter *T)
{
    // Clones share table/index with the parent — only free the wrapper.
    // The parent (allocated on the caller's stack) is freed via
    // free_kmer_adapter, which releases table/index.
    if (T->clone)
        free(T);
}

uint8 *adapter_Current_Entry(KmerStreamAdapter *T, uint8 *ent)
{
    if (ent == NULL) {
        ent = (uint8 *) malloc(T->tbyte);
        if (ent == NULL) { perror("malloc"); exit(1); }
        if (T->csuf == NULL) return ent;
    }
    // ibyte=3 prefix bytes (big-endian, matching libfastk.c case 3)
    ent[0] = (T->cpre >> 16) & 0xff;
    ent[1] = (T->cpre >>  8) & 0xff;
    ent[2] =  T->cpre        & 0xff;
    memcpy(ent + T->ibyte, T->csuf, T->pbyte);
    return ent;
}