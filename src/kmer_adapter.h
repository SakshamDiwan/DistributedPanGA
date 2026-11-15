#ifndef KMER_ADAPTER_H
#define KMER_ADAPTER_H

#include <stdint.h>
#include "extract.h"    // for RecordBuffer
#include "record.h"     // for RecordSizing, offsets
#include "constants.h"  // for KBYTES

// KmerStreamAdapter — our sorted RecordBuffer as a FastK-style k-mer table,
// walked by the merge stage through the adapter_* functions below. Field names
// follow FastK's Kmer_Stream, but this is our own type: nothing casts it to
// Kmer_Stream, so field order and widths are free to change.

typedef struct {
    int     kmer;       // k-mer length in bases (40)
    int64   nels;       // number of table entries (= records in this rank's slice)
    int64   cidx;       // cursor: index of the current entry
    uint8  *csuf;       // cursor: pointer to the current entry in table
    int     cpre;       // cursor: 24-bit prefix (first ibyte k-mer bytes) of the current entry
    int     ibyte;      // k-mer bytes turned into the prefix, not stored per entry (3)
    int     kbyte;      // bytes of the whole k-mer (10)
    int     tbyte;      // prefix + entry bytes, as Current_Entry returns it (ibyte + pbyte)
    int     hbyte;      // k-mer bytes stored per entry (kbyte - ibyte)
    int     pbyte;      // bytes per table entry: hbyte + mask + LCP + payload
    int     ixlen;      // number of possible prefixes = length of index (2^(8*ibyte))
    uint8  *table;      // nels entries of pbyte bytes, sorted by k-mer
    int64  *index;      // index[p] = number of entries with prefix <= p
    uint8  *ctop;       // one past the last entry of table
    int     clone;      // 1 if a per-thread copy sharing table/index (free only the wrapper)
} KmerStreamAdapter;

// PostListAdapter — the payload description FastGA's merge needs (field names
// follow FastGA's Post_List; our own type, never cast).

typedef struct {
    int     pbyte;         // payload bytes per entry: position + contig/strand
    int     cbyte;         // contig/strand bytes (strand = top bit of the last byte)
    int64   maxp;          // largest prefix group + 1 (sizes the merge thread's cache)
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