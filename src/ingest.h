#ifndef INGEST_H
#define INGEST_H

#include <stdint.h>
#include "GDB.h"

// Record layout for one k-mer entry:
//   [LCP: 1 byte] [packed k-mer: KBYTES bytes] [mask: 1 byte]
//   [position: PostBytes bytes] [contig+strand: ContBytes bytes]
//
// LCP is at byte 0. msd_sort starts at digit=1 (byte 1 = first kmer byte),
// so it sorts the kmer correctly and writes inter-bucket LCP values to byte 0.
// We discard msd_sort's LCP writes and recompute LCP ourselves post-sort.
// For k=40: KBYTES=10, so record is 1 + 10 + 1 + PostBytes + ContBytes

typedef struct {
    uint8_t *records;     // buffer of packed k-mer records
    int64_t  count;       // number of records in buffer
    int64_t  capacity;    // allocated capacity
    int      record_size; // bytes per record
} KmerBuffer;

// Initialize a k-mer buffer
void kmer_buffer_init(KmerBuffer *buf, int record_size, int64_t initial_capacity);

// Free a k-mer buffer
void kmer_buffer_free(KmerBuffer *buf);

// Extract syncmer-filtered k-mers from assigned contigs
//   gdb        - opened genome database
//   rank       - this MPI rank
//   num_ranks  - total MPI ranks
//   kmer       - k-mer size (40)
//   buf        - output buffer (filled with k-mer records)
//   counts     - caller-allocated int64_t[1024], filled with per-bucket k-mer counts
void extract_kmers(GDB *gdb, int rank, int num_ranks, int kmer, KmerBuffer *buf, int64_t *counts);

#endif