#ifndef KMER_ADAPTER_H
#define KMER_ADAPTER_H

#include <stdint.h>
#include "extract.h"    // for RecordBuffer
#include "record.h"     // for RecordSizing, offsets
#include "constants.h"  // for KBYTES

// KmerStreamAdapter — wraps our RecordBuffer to look like Kmer_Stream.

typedef struct {
    // ---- Mirror of Kmer_Stream (libfastk.h lines 65..93) ----
    int     kmer;
    int     minval;
    int64   nels;
    int64   cidx;
    uint8  *csuf;
    int     cpre;
    int     ibyte;
    int     kbyte;
    int     tbyte;
    int     hbyte;
    int     pbyte;
    int     ixlen;
    int     shift;
    uint8  *table;
    int64  *index;
    int    *inver;
    int     copn;       // unused by adapter, kept for layout compat
    int     part;       // unused
    int     nthr;       // unused
    int     nlen;       // unused
    char   *name;       // unused
    uint8  *ctop;
    int64  *neps;       // unused
    int     clone;
    // ---- end mirror ----
} KmerStreamAdapter;

// PostListAdapter — wraps post info to look like Post_List.
// Fields must match Post_List layout (first N fields in order).

typedef struct {
    int     pbyte;         // = PostBytes + ContBytes
    int     cbyte;         // = ContBytes
    int64   nels;          // total number of positions
    int64   maxp;          // max posts per k-mer prefix (for cache sizing)
    int     freq;          // frequency cutoff (not used in self-merge)
    int     nctg;          // number of contigs
    int    *perm;          // contig permutation
    // We don't need the rest (index, cache, file handles, etc.)
} PostListAdapter;

int  build_kmer_adapter(const RecordBuffer *buf,
                        KmerStreamAdapter *T,
                        PostListAdapter *P);

void free_kmer_adapter(KmerStreamAdapter *T, PostListAdapter *P);

void adapter_First_Kmer_Entry(KmerStreamAdapter *T);
void adapter_Next_Kmer_Entry(KmerStreamAdapter *T);
void adapter_GoTo_Kmer_Index(KmerStreamAdapter *T, int64 i);
KmerStreamAdapter *adapter_Clone_Kmer_Stream(KmerStreamAdapter *T);
void adapter_Free_Kmer_Stream(KmerStreamAdapter *T);
uint8 *adapter_Current_Entry(KmerStreamAdapter *T, uint8 *ent);

#endif