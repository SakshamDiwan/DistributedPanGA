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

// FastK names, for reading upstream code: kmer, nels, cidx, csuf, cpre, ibyte,
// kbyte, tbyte, hbyte, pbyte, ixlen, table, index, ctop, clone (in field order).
typedef struct {
    int     kmer_length_bases;          // k-mer length in bases (40)
    int64   num_entries;                // table entries (= records in this rank's slice)
    int64   current_entry_index;        // cursor: index of the current entry
    uint8  *current_entry;              // cursor: pointer to the current entry in entries
    int     current_prefix;             // cursor: 24-bit prefix (group) of the current entry
    int     prefix_bytes;               // leading k-mer bytes turned into the prefix, not stored per entry (3)
    int     kmer_bytes;                 // bytes of a whole k-mer, 2 bits per base (10)
    int     prefix_plus_entry_bytes;    // prefix + entry, as Current_Entry rebuilds it (16)
    int     entry_kmer_suffix_bytes;    // k-mer bytes stored per entry: kmer_bytes - prefix_bytes (7)
    int     entry_bytes;                // bytes per entry: k-mer suffix + mask + LCP + payload (13)
    int     num_prefixes;               // possible prefixes = length of the cumulative count (2^24)
    uint8  *entries;                    // num_entries entries of entry_bytes bytes, sorted by k-mer
    int64  *cumulative_count_by_prefix; // [p] = number of entries with prefix <= p
    uint8  *entries_end;                // one past the last entry
    int     is_clone;                   // per-thread copy sharing entries/counts (free only the wrapper)
} KmerStreamAdapter;

// PostListAdapter — the payload description FastGA's merge needs (field names
// follow FastGA's Post_List; our own type, never cast).

// FastGA names: pbyte, cbyte, maxp (in field order).
typedef struct {
    int     payload_bytes;         // payload bytes per entry: position + contig/strand (4)
    int     contig_strand_bytes;   // contig/strand bytes; strand = top bit of the last byte (1)
    int64   largest_prefix_group_num_entries_plus_one;
                                   // entries in the largest prefix group, + 1. Sizes the
                                   // merge thread's cache, which adds its own +1 for the
                                   // end-of-group marker -- so this +1 is redundant
                                   // (follow-up: store the plain count, update test_adapter).
} PostListAdapter;

int  build_kmer_adapter(const RecordBuffer *buffer,
                        KmerStreamAdapter *kmer_table,
                        PostListAdapter *payload_layout);

void free_kmer_adapter(KmerStreamAdapter *kmer_table, PostListAdapter *payload_layout);

void adapter_First_Kmer_Entry(KmerStreamAdapter *kmer_table);
void adapter_Next_Kmer_Entry(KmerStreamAdapter *kmer_table);
void adapter_GoTo_Kmer_Index(KmerStreamAdapter *kmer_table, int64 i);
KmerStreamAdapter *adapter_Clone_Kmer_Stream(KmerStreamAdapter *kmer_table);
void adapter_Free_Kmer_Stream(KmerStreamAdapter *kmer_table);
uint8 *adapter_Current_Entry(KmerStreamAdapter *kmer_table, uint8 *ent);

#endif