#ifndef DPGA_EXTRACT_H
#define DPGA_EXTRACT_H

#include <stdint.h>
#include "GDB.h"
#include "record.h"

// --- Phase 1 path: TupleBuffer (struct of fields, easy to inspect) -----------

typedef struct {
    int     contig;
    int64_t position;     // j (forward) or j+TMER (RC)
    int     strand;       // 0 = forward, 1 = RC
    uint8_t kmer[10];     // KBYTES, MSB-first 2-bit packing
} Tuple;

typedef struct {
    Tuple   *data;
    int64_t  count;
    int64_t  capacity;
} TupleBuffer;

void tuple_buffer_init(TupleBuffer *b);
void tuple_buffer_free(TupleBuffer *b);

void extract_to_tuples(GDB *gdb, int num_workers, int worker_id, TupleBuffer *buf);

// --- Phase 2 path: RecordBuffer (flat byte records, msd_sort-ready) ----------

typedef struct {
    uint8_t      *data;       // contiguous record bytes; NOTE: trailing sentinel byte
    int64_t       count;      // number of records currently in buffer
    int64_t       capacity;   // capacity in records (not bytes)
    RecordSizing  sizing;
} RecordBuffer;

void record_buffer_init(RecordBuffer *b, RecordSizing sizing);
void record_buffer_free(RecordBuffer *b);

void extract_to_records(GDB *gdb, int num_workers, int worker_id, RecordBuffer *buf);

#endif
