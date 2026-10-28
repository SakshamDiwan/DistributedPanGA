// ============================================================================
// R-06 (MPI tier) -- chunked_alltoallv end to end.  See tests/CATALOG.md.
// ============================================================================
//
// test_chunked_exchange checks the round *plan* without MPI. This test runs
// the real chunked_alltoallv between ranks and compares the received bytes with
// one plain MPI_Alltoallv of the same data. It covers:
//   - memory buffers (stage-2 style), memory callbacks, and file callbacks with
//     pread/pwrite (the seed-streaming path in redistribute_seeds);
//   - record sizes 16 and 9;
//   - budgets from one record per peer (many rounds) to large (direct path);
//   - a DIFFERENT budget on every rank, which used to risk mismatched
//     collectives: the ranks must agree on the minimum.
// Data sizes are uneven: some peers get nothing, one pair gets much more.
//
// Run under srun with any world size; exit 0 iff every check passed on every
// rank.

#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>

#include "chunked_exchange.h"

static int rank, nranks, failures = 0, checks = 0;

static uint8_t byte_at(int src, int dst, int64_t off)
{
    uint64_t x = (uint64_t) off * 0x9E3779B97F4A7C15ULL ^ ((uint64_t) src << 40)
               ^ ((uint64_t) dst << 20);
    x ^= x >> 29;
    return (uint8_t) (x * 0xBF58476D1CE4E5B9ULL >> 56);
}

static int64_t records_for(int src, int dst)
{
    if (src != dst && (src * 7 + dst * 3) % 5 == 0) return 0;   // some empty pairs
    if (src == 0 && dst == nranks - 1) return 900;       // one heavy pair
    return 1 + (src * 31 + dst * 17) % 97;
}

// --- callback sources/sinks -----------------------------------------------

typedef struct { uint8_t **per_peer; } MemSrc;
static void mem_fill(void *c, int p, int64_t off, int64_t len, uint8_t *dst)
{ memcpy(dst, ((MemSrc *) c)->per_peer[p] + off, len); }

typedef struct { uint8_t *buf; const int64_t *prefix; } MemSink;
static void mem_drain(void *c, int p, int64_t off, int64_t len, const uint8_t *src)
{ MemSink *s = c; memcpy(s->buf + s->prefix[p] + off, src, len); }

typedef struct { int *fd; } FileSrc;
static void file_fill(void *c, int p, int64_t off, int64_t len, uint8_t *dst)
{
    int64_t got = 0;
    while (got < len) {
        ssize_t n = pread(((FileSrc *) c)->fd[p], dst + got, len - got, off + got);
        if (n <= 0) { perror("pread"); MPI_Abort(MPI_COMM_WORLD, 1); }
        got += n;
    }
}

typedef struct { int fd; const int64_t *prefix; } FileSink;
static void file_drain(void *c, int p, int64_t off, int64_t len, const uint8_t *src)
{
    FileSink *s = c;
    int64_t put = 0;
    while (put < len) {
        ssize_t n = pwrite(s->fd, src + put, len - put, s->prefix[p] + off + put);
        if (n <= 0) { perror("pwrite"); MPI_Abort(MPI_COMM_WORLD, 1); }
        put += n;
    }
}

// --- one configuration ----------------------------------------------------

enum { MODE_BUF, MODE_MEMCB, MODE_FILECB };

static void run_case(int unit, int64_t budget_mine, int mode, int want_multi,
                     const char *label)
{
    int64_t *sbytes = calloc(nranks, sizeof(int64_t));
    int64_t *rbytes = calloc(nranks, sizeof(int64_t));
    int64_t *sprefix = calloc(nranks, sizeof(int64_t));
    int64_t *rprefix = calloc(nranks, sizeof(int64_t));
    int64_t stot = 0, rtot = 0;

    for (int p = 0; p < nranks; p++) {
        sbytes[p] = records_for(rank, p) * unit;
        sprefix[p] = stot;
        stot += sbytes[p];
    }
    uint8_t *send = malloc(stot + 1);
    uint8_t **per_peer = malloc(nranks * sizeof(uint8_t *));
    for (int p = 0; p < nranks; p++) {
        per_peer[p] = send + sprefix[p];
        for (int64_t i = 0; i < sbytes[p]; i++) per_peer[p][i] = byte_at(rank, p, i);
    }

    xchg_exchange_sizes(sbytes, rbytes, MPI_COMM_WORLD);
    for (int p = 0; p < nranks; p++) {
        rprefix[p] = rtot;
        rtot += rbytes[p];
    }

    // Reference: one plain Alltoallv (sizes here are tiny, so ints are fine).
    int *sc = malloc(nranks * sizeof(int)), *sd = malloc(nranks * sizeof(int));
    int *rc = malloc(nranks * sizeof(int)), *rd = malloc(nranks * sizeof(int));
    for (int p = 0; p < nranks; p++) {
        sc[p] = (int) sbytes[p]; sd[p] = (int) sprefix[p];
        rc[p] = (int) rbytes[p]; rd[p] = (int) rprefix[p];
    }
    uint8_t *ref = malloc(rtot + 1), *got = calloc(rtot + 1, 1);
    MPI_Alltoallv(send, sc, sd, MPI_BYTE, ref, rc, rd, MPI_BYTE, MPI_COMM_WORLD);

    XchgStats st;
    char tmpl[64];
    int out_fd = -1, *in_fd = NULL;
    MemSrc msrc = { per_peer };
    MemSink mdst = { got, rprefix };
    XchgSource src = { NULL, NULL, NULL };
    XchgSink dst = { NULL, NULL, NULL };
    FileSrc fsrc;
    FileSink fdst;

    if (mode == MODE_BUF) {
        src.buf = send;
        dst.buf = got;
    } else if (mode == MODE_MEMCB) {
        src.fill = mem_fill;  src.ctx = &msrc;
        dst.drain = mem_drain; dst.ctx = &mdst;
    } else {
        in_fd = malloc(nranks * sizeof(int));
        for (int p = 0; p < nranks; p++) {
            snprintf(tmpl, sizeof(tmpl), "/tmp/xchg-in-%d-%d-XXXXXX", rank, p);
            in_fd[p] = mkstemp(tmpl);
            unlink(tmpl);
            if (sbytes[p] && pwrite(in_fd[p], per_peer[p], sbytes[p], 0) != sbytes[p]) {
                perror("pwrite in"); MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        snprintf(tmpl, sizeof(tmpl), "/tmp/xchg-out-%d-XXXXXX", rank);
        out_fd = mkstemp(tmpl);
        unlink(tmpl);
        fsrc.fd = in_fd;
        fdst.fd = out_fd; fdst.prefix = rprefix;
        src.fill = file_fill;  src.ctx = &fsrc;
        dst.drain = file_drain; dst.ctx = &fdst;
    }

    chunked_alltoallv(sbytes, rbytes, unit, budget_mine, &src, &dst, MPI_COMM_WORLD, &st);

    if (mode == MODE_FILECB && rtot > 0 && pread(out_fd, got, rtot, 0) != rtot) {
        perror("pread out"); MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int ok = (memcmp(ref, got, rtot) == 0);
    int64_t min_budget;
    MPI_Allreduce(&budget_mine, &min_budget, 1, MPI_INT64_T, MPI_MIN, MPI_COMM_WORLD);
    int64_t want_cap = xchg_chunk_cap(nranks, unit, min_budget);
    if (!st.direct && st.cap != want_cap) ok = 0;    // all ranks used the minimum budget
    if (want_multi && (st.direct || st.rounds < 2)) ok = 0;   // must really chunk

    int all_ok;
    MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    checks++;
    if (!all_ok) failures++;
    if (rank == 0)
        printf("  %-44s rounds=%-5lld %-7s %s\n", label, (long long) st.rounds,
               st.direct ? "direct" : "chunked", all_ok ? "ok" : "FAIL");

    if (in_fd) { for (int p = 0; p < nranks; p++) close(in_fd[p]); free(in_fd); }
    if (out_fd >= 0) close(out_fd);
    free(sbytes); free(rbytes); free(sprefix); free(rprefix);
    free(send); free(per_peer); free(ref); free(got);
    free(sc); free(sd); free(rc); free(rd);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
    if (rank == 0) printf("test_chunked_mpi (world size %d)\n", nranks);

    const char *mode_name[] = { "buffers", "memory callbacks", "file pread/pwrite" };
    for (int mode = 0; mode < 3; mode++) {
        for (int u = 0; u < 2; u++) {
            int unit = u ? 9 : 16;
            char label[96];
            // one record per peer per round on the smallest-budget rank,
            // and a different (larger) budget on every other rank
            snprintf(label, sizeof(label), "%s, unit %d, tiny mixed budgets",
                     mode_name[mode], unit);
            run_case(unit, (int64_t) unit * nranks * (rank + 1), mode, 1, label);
            snprintf(label, sizeof(label), "%s, unit %d, 4 KiB budget", mode_name[mode], unit);
            run_case(unit, 4096, mode, 1, label);
            snprintf(label, sizeof(label), "%s, unit %d, default budget", mode_name[mode], unit);
            run_case(unit, XCHG_DEFAULT_ROUND_BYTES, mode, 0, label);
        }
    }
    if (rank == 0)
        printf("test_chunked_mpi: %d checks, %d failures\n", checks, failures);
    MPI_Finalize();
    return failures ? 1 : 0;
}
