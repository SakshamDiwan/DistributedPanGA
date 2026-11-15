#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kmer_adapter.h"

int build_kmer_adapter(const RecordBuffer *buffer, KmerStreamAdapter *kmer_table, PostListAdapter *payload_layout)
{
    int64 num_records = buffer->count;
    int record_size = buffer->sizing.record_size;
    int position_bytes = buffer->sizing.position_bytes;
    int contig_strand_bytes = buffer->sizing.contig_strand_bytes;

    kmer_table->kmer_length_bases       = 40;
    kmer_table->prefix_bytes            = 3;
    kmer_table->kmer_bytes              = 10;
    kmer_table->entry_kmer_suffix_bytes = kmer_table->kmer_bytes - kmer_table->prefix_bytes;   // 7
    kmer_table->entry_bytes             = kmer_table->entry_kmer_suffix_bytes + 1 + 1          // + mask + LCP
                                          + position_bytes + contig_strand_bytes;             // 13
    kmer_table->prefix_plus_entry_bytes = kmer_table->prefix_bytes + kmer_table->entry_bytes; // 16, used by Current_Entry
    kmer_table->num_entries             = num_records;
    kmer_table->num_prefixes            = 1 << (8 * kmer_table->prefix_bytes);               // 2^24

    kmer_table->entries = (uint8 *) malloc(num_records * kmer_table->entry_bytes);
    if (kmer_table->entries == NULL) {
        fprintf(stderr, "build_kmer_adapter: out of memory for table\n");
        return -1;
    }

    kmer_table->cumulative_count_by_prefix = (int64 *) calloc(kmer_table->num_prefixes, sizeof(int64));
    if (kmer_table->cumulative_count_by_prefix == NULL) {
        fprintf(stderr, "build_kmer_adapter: out of memory for index\n");
        free(kmer_table->entries);
        return -1;
    }

    uint8 *dst = kmer_table->entries;
    int64 max_run = 0, current_run = 0;
    int   prev_prefix = -1;

    for (int64 i = 0; i < num_records; i++) {
        const uint8 *src = buffer->data + i * record_size;

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

        memcpy(dst, src + KMER_OFFSET + kmer_table->prefix_bytes, kmer_table->entry_kmer_suffix_bytes);
        dst += kmer_table->entry_kmer_suffix_bytes;
        *dst++ = 0;
        *dst++ = src[LCP_OFFSET];
        memcpy(dst, src + POS_OFFSET, position_bytes + contig_strand_bytes);
        dst += position_bytes + contig_strand_bytes;
    }
    if (current_run > max_run) max_run = current_run;

    // "≤ p" cumulative: cumulative_count_by_prefix[p] = number of records with prefix ≤ p.
    // Matches libfastk.c semantics (see First_Kmer_Entry, Next_Kmer_Entry).
    memset(kmer_table->cumulative_count_by_prefix, 0, kmer_table->num_prefixes * sizeof(int64));
    for (int64 i = 0; i < num_records; i++) {
        const uint8 *src = buffer->data + i * record_size;
        int prefix = (src[KMER_OFFSET] << 16) |
                    (src[KMER_OFFSET + 1] << 8) |
                    src[KMER_OFFSET + 2];
        kmer_table->cumulative_count_by_prefix[prefix]++;
    }
    {
        int64 cum = 0;
        for (int64 p = 0; p < kmer_table->num_prefixes; p++) {
            cum += kmer_table->cumulative_count_by_prefix[p];
            kmer_table->cumulative_count_by_prefix[p] = cum;
        }
    }

    kmer_table->current_entry       = kmer_table->entries;
    kmer_table->current_entry_index = 0;
    kmer_table->current_prefix      = 0;
    while (kmer_table->current_prefix < kmer_table->num_prefixes
           && kmer_table->cumulative_count_by_prefix[kmer_table->current_prefix] <= 0)
        kmer_table->current_prefix++;
    kmer_table->is_clone    = 0;
    kmer_table->entries_end = kmer_table->entries + num_records * kmer_table->entry_bytes;

    payload_layout->payload_bytes       = position_bytes + contig_strand_bytes;
    payload_layout->contig_strand_bytes = contig_strand_bytes;
    payload_layout->largest_prefix_group_num_entries_plus_one = max_run + 1;

    return 0;
}

void free_kmer_adapter(KmerStreamAdapter *kmer_table, PostListAdapter *payload_layout)
{
    (void) payload_layout;
    if (kmer_table->entries) { free(kmer_table->entries); kmer_table->entries = NULL; }
    if (kmer_table->cumulative_count_by_prefix) { free(kmer_table->cumulative_count_by_prefix); kmer_table->cumulative_count_by_prefix = NULL; }
}

void adapter_First_Kmer_Entry(KmerStreamAdapter *kmer_table)
{
    kmer_table->current_entry = kmer_table->entries;
    kmer_table->current_entry_index = 0;
    kmer_table->current_prefix = 0;
    while (kmer_table->current_prefix < kmer_table->num_prefixes && kmer_table->cumulative_count_by_prefix[kmer_table->current_prefix] <= 0)
        kmer_table->current_prefix++;
}

void adapter_Next_Kmer_Entry(KmerStreamAdapter *kmer_table)
{
    kmer_table->current_entry += kmer_table->entry_bytes;
    kmer_table->current_entry_index++;
    if (kmer_table->current_entry >= kmer_table->entries_end) {
        kmer_table->current_entry = NULL;
        kmer_table->current_prefix = kmer_table->num_prefixes;
        return;
    }
    while (kmer_table->current_prefix < kmer_table->num_prefixes && kmer_table->cumulative_count_by_prefix[kmer_table->current_prefix] <= kmer_table->current_entry_index)
        kmer_table->current_prefix++;
}

// O(num_prefixes) linear scan — acceptable only because self_adaptamer_merge
// only calls this NTHREADS-1 times. If NTHREADS grows, build kmer_table->inver
// in build_kmer_adapter and use S->inver[i>>S->shift] like libfastk.c:1289.
void adapter_GoTo_Kmer_Index(KmerStreamAdapter *kmer_table, int64 i)
{
    if (i >= kmer_table->num_entries) {
        kmer_table->current_entry = NULL;
        kmer_table->current_prefix = kmer_table->num_prefixes;
        kmer_table->current_entry_index = kmer_table->num_entries;
        return;
    }
    kmer_table->current_entry_index = i;
    kmer_table->current_entry = kmer_table->entries + i * kmer_table->entry_bytes;
    kmer_table->current_prefix = 0;
    while (kmer_table->current_prefix < kmer_table->num_prefixes && kmer_table->cumulative_count_by_prefix[kmer_table->current_prefix] <= i)
        kmer_table->current_prefix++;
}

KmerStreamAdapter *adapter_Clone_Kmer_Stream(KmerStreamAdapter *kmer_table)
{
    KmerStreamAdapter *copy = malloc(sizeof(KmerStreamAdapter));
    *copy = *kmer_table;
    copy->is_clone = 1;
    return copy;
}

void adapter_Free_Kmer_Stream(KmerStreamAdapter *kmer_table)
{
    // Clones share entries/cumulative_count_by_prefix with the parent — only free the wrapper.
    // The parent (allocated on the caller's stack) is freed via
    // free_kmer_adapter, which releases entries/cumulative_count_by_prefix.
    if (kmer_table->is_clone)
        free(kmer_table);
}

uint8 *adapter_Current_Entry(KmerStreamAdapter *kmer_table, uint8 *ent)
{
    if (ent == NULL) {
        ent = (uint8 *) malloc(kmer_table->prefix_plus_entry_bytes);
        if (ent == NULL) { perror("malloc"); exit(1); }
        if (kmer_table->current_entry == NULL) return ent;
    }
    // prefix_bytes=3 prefix bytes (big-endian, matching libfastk.c case 3)
    ent[0] = (kmer_table->current_prefix >> 16) & 0xff;
    ent[1] = (kmer_table->current_prefix >>  8) & 0xff;
    ent[2] =  kmer_table->current_prefix        & 0xff;
    memcpy(ent + kmer_table->prefix_bytes, kmer_table->current_entry, kmer_table->entry_bytes);
    return ent;
}