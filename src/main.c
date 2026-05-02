#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "GDB.h"
#include "ingest.h"

extern void msd_sort(uint8_t *array, int64_t nelem, int rsize, int ksize,
                     int64_t *part, int beg, int end, int nthreads);

int main(int argc, char *argv[])
{
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2)
    {   if (rank == 0)
            fprintf(stderr, "Usage: distpanga <genome.1gdb>\n");
        MPI_Finalize();
        return 1;
    }

    // All ranks read the GDB (metadata is small)
    GDB _gdb, *gdb = &_gdb;
    Read_GDB(gdb, argv[1]);

    if (rank == 0)
        printf("Loaded GDB: %d contigs, %lld total bases\n", gdb->ncontig, gdb->seqtot);

    // Compute record sizing (must match ingest.c's extract_kmers exactly)
    int kbytes     = 40 / 4;  // GIXmake stores 9 kmer bytes (neq[bost+4..bost+36]), not 10
    int post_bytes = 0;
    { int64_t cum = 1;
      while (cum < gdb->maxctg)
      { cum *= 256; post_bytes++; }
    }
    int cont_bytes = 0;
    { int64_t cum = 1;
      while (cum < 2 * gdb->ncontig)
      { cum *= 256; cont_bytes++; }
    }
    int record_size = 1 + kbytes + 1 + post_bytes + cont_bytes; // LCP + kmer + mask + pos + cont

    // Extract k-mers
    KmerBuffer buf;
    int64_t local_counts[1024];
    extract_kmers(gdb, rank, size, 40, &buf, local_counts);

    // Sum up total k-mers across all ranks (for reporting)
    int64_t total_kmers;
    MPI_Reduce(&buf.count, &total_kmers, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
        printf("Total k-mers extracted: %lld\n", total_kmers);

    // Every rank gets the global histogram (all ranks need it to compute Ksplit)
    int64_t global_counts[1024];
    MPI_Allreduce(local_counts, global_counts, 1024, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    // Print bucket counts for comparison against GIXmake DEBUG_SETUP output
    if (rank == 0)
    {   printf("\nBuckets\n");
        for (int i = 0; i < 1024; i++)
            printf(" %3d: %10lld\n", i, global_counts[i]);
        fflush(stdout);
    }

    if (rank == 0)
    {   FILE *df = fopen("ours_bucket870.txt", "w");
        if (df != NULL)
        {   int64_t n_dumped = 0;
            for (int64_t i = 0; i < buf.count; i++)
            {   uint8_t *rec = buf.records + i * record_size;
                int bucket = (rec[1] << 2) | (rec[2] >> 6);
                if (bucket == 870) n_dumped++;
            }
            fprintf(df, "# bucket 870: %lld records, record_size=%d\n",
                    n_dumped, record_size);
            for (int64_t i = 0; i < buf.count; i++)
            {   uint8_t *rec = buf.records + i * record_size;
                int bucket = (rec[1] << 2) | (rec[2] >> 6);
                if (bucket != 870) continue;
                for (int q = 0; q < record_size; q++)
                    fprintf(df, "%02x", rec[q]);
                fprintf(df, "\n");
            }
            fclose(df);
        }
    }

    // Convert global_counts to prefix sums (same as GIXmake)
    for (int i = 1; i < 1024; i++)
        global_counts[i] += global_counts[i-1];

    // Compute Ksplit[0..size] and Select[0..1023] — identical on every rank
    int *Ksplit = (int *) malloc((size + 1) * sizeof(int));
    int  Select[1024];

    Ksplit[0] = 0;
    int    n = 1;
    int64_t t = global_counts[1023] / size;
    for (int i = 0; i < 1024; i++)
    {   if (global_counts[i] >= t)
        {   if (i > 0 && global_counts[i] - t > t - global_counts[i-1])
            {   Select[i] = n;
                Ksplit[n] = i;
            }
            else
            {   Select[i] = n - 1;
                Ksplit[n] = i + 1;
            }
            n += 1;
            t = ((int64_t)n * global_counts[1023]) / size;
        }
        else
            Select[i] = n - 1;
    }
    Ksplit[size] = 1024;

    // Step 1: count how many records this rank sends to each destination rank
    int *send_counts = (int *) calloc(size, sizeof(int));
    for (int64_t i = 0; i < buf.count; i++)
    {   uint8_t *rec = buf.records + i * record_size;
        int bucket = (rec[1] << 2) | (rec[2] >> 6);
        send_counts[Select[bucket]]++;
    }

    // Convert record counts to byte counts for MPI
    int *send_bytes = (int *) malloc(size * sizeof(int));
    for (int i = 0; i < size; i++)
        send_bytes[i] = send_counts[i] * record_size;

    // Step 2: build contiguous send buffer grouped by destination rank
    int *send_displs = (int *) malloc(size * sizeof(int));
    send_displs[0] = 0;
    for (int i = 1; i < size; i++)
        send_displs[i] = send_displs[i-1] + send_bytes[i-1];

    uint8_t *send_buf = (uint8_t *) malloc(buf.count * record_size);
    int *cursors = (int *) calloc(size, sizeof(int));  // current write offset per rank
    for (int64_t i = 0; i < buf.count; i++)
    {   uint8_t *rec = buf.records + i * record_size;
        int bucket = (rec[1] << 2) | (rec[2] >> 6);
        int dest   = Select[bucket];
        uint8_t *dst = send_buf + send_displs[dest] + cursors[dest];
        memcpy(dst, rec, record_size);
        cursors[dest] += record_size;
    }
    free(cursors);
    kmer_buffer_free(&buf);

    // Step 3: exchange byte counts so every rank knows how much it will receive
    int *recv_bytes  = (int *) malloc(size * sizeof(int));
    MPI_Alltoall(send_bytes, 1, MPI_INT, recv_bytes, 1, MPI_INT, MPI_COMM_WORLD);

    int *recv_displs = (int *) malloc(size * sizeof(int));
    recv_displs[0] = 0;
    for (int i = 1; i < size; i++)
        recv_displs[i] = recv_displs[i-1] + recv_bytes[i-1];

    int64_t total_recv_bytes = recv_displs[size-1] + recv_bytes[size-1];
    uint8_t *recv_buf = (uint8_t *) malloc(total_recv_bytes + 1); // +1 for msd_sort sentinel

    // Step 4: redistribute k-mers to their sort-owner ranks
    MPI_Alltoallv(send_buf,  send_bytes,  send_displs, MPI_BYTE,
                  recv_buf,  recv_bytes,  recv_displs, MPI_BYTE,
                  MPI_COMM_WORLD);

    int64_t recv_count = total_recv_bytes / record_size;
    printf("Rank %d: received %lld k-mers for sorting\n", rank, recv_count);
    fflush(stdout);

    free(send_buf);
    free(send_bytes);
    free(send_counts);
    free(send_displs);
    free(recv_bytes);
    free(recv_displs);
    free(recv_buf);
    free(Ksplit);

    Close_GDB(gdb);
    MPI_Finalize();
    return 0;
}