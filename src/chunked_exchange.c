#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chunked_exchange.h"

int64_t xchg_round_budget_from_env(void)
{
    const char *s = getenv("PGA_EXCHANGE_ROUND_BYTES");
    char *end;
    long long v;

    if (s == NULL || *s == '\0')
        return XCHG_DEFAULT_ROUND_BYTES;
    v = strtoll(s, &end, 10);
    if (end == s || v <= 0)
        return XCHG_DEFAULT_ROUND_BYTES;
    return (int64_t) v;
}

int64_t xchg_chunk_cap(int npeers, int unit, int64_t budget)
{
    int64_t cap, hard;

    if (npeers < 1) npeers = 1;
    if (unit < 1) unit = 1;
    cap  = (budget / npeers) / unit * unit;
    hard = ((int64_t) INT_MAX / npeers) / unit * unit;
    if (cap > hard) cap = hard;
    if (cap < unit) cap = unit;
    return cap;
}

int64_t xchg_rounds_for(int64_t bytes, int64_t cap)
{
    if (bytes <= 0) return 0;
    return (bytes + cap - 1) / cap;
}

int64_t xchg_round_count(int64_t bytes, int64_t cap, int64_t k)
{
    int64_t off = k * cap;
    if (off >= bytes) return 0;
    return (bytes - off < cap) ? bytes - off : cap;
}

int64_t xchg_plan_round(const int64_t *bytes, int npeers, int unit, int64_t cap,
                        int64_t k, int *counts, int *displs)
{
    int64_t total = 0;
    for (int p = 0; p < npeers; p++)
    {
        if (bytes[p] < 0) return -1;
        if (bytes[p] % unit != 0) return -2;
        int64_t c = xchg_round_count(bytes[p], cap, k);
        if (total + c > INT_MAX) return -3;
        counts[p] = (int) c;
        displs[p] = (int) total;
        total += c;
    }
    return total;
}

void xchg_exchange_sizes(const int64_t *send_bytes, int64_t *recv_bytes, MPI_Comm comm)
{
    MPI_Alltoall((void *) send_bytes, 1, MPI_INT64_T, recv_bytes, 1, MPI_INT64_T, comm);
}

static void die(MPI_Comm comm, const char *what, long long n)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    fprintf(stderr, "chunked_alltoallv: rank %d: %s (%lld)\n", rank, what, n);
    MPI_Abort(comm, 1);
}

void chunked_alltoallv(const int64_t *send_bytes, const int64_t *recv_bytes,
                       int unit, int64_t budget,
                       const XchgSource *src, const XchgSink *dst,
                       MPI_Comm comm, XchgStats *stats)
{
    int n;
    MPI_Comm_size(comm, &n);

    int64_t *sprefix = malloc(n * sizeof(int64_t));
    int64_t *rprefix = malloc(n * sizeof(int64_t));
    int *scounts = malloc(n * sizeof(int)), *sdispls = malloc(n * sizeof(int));
    int *rcounts = malloc(n * sizeof(int)), *rdispls = malloc(n * sizeof(int));
    if (!sprefix || !rprefix || !scounts || !sdispls || !rcounts || !rdispls)
        die(comm, "out of memory for count arrays", n);

    int64_t stot = 0, rtot = 0;
    for (int p = 0; p < n; p++)
    {
        if (send_bytes[p] < 0 || send_bytes[p] % unit) die(comm, "bad send size", send_bytes[p]);
        if (recv_bytes[p] < 0 || recv_bytes[p] % unit) die(comm, "bad receive size", recv_bytes[p]);
        sprefix[p] = stot; stot += send_bytes[p];
        rprefix[p] = rtot; rtot += recv_bytes[p];
    }

    XchgStats st;
    memset(&st, 0, sizeof(st));
    st.sent_bytes = stot;
    st.recv_bytes = rtot;

    // Direct path: one Alltoallv, exactly as before, when it is safe everywhere.
    int64_t mine = (stot > rtot) ? stot : rtot, worst;
    MPI_Allreduce(&mine, &worst, 1, MPI_INT64_T, MPI_MAX, comm);
    int contiguous = (src->buf != NULL && dst->buf != NULL), all_contig;
    MPI_Allreduce(&contiguous, &all_contig, 1, MPI_INT, MPI_MIN, comm);
    if (all_contig && worst <= budget && worst <= INT_MAX)
    {
        for (int p = 0; p < n; p++)
        {
            scounts[p] = (int) send_bytes[p]; sdispls[p] = (int) sprefix[p];
            rcounts[p] = (int) recv_bytes[p]; rdispls[p] = (int) rprefix[p];
        }
        MPI_Alltoallv((void *) src->buf, scounts, sdispls, MPI_BYTE,
                      dst->buf, rcounts, rdispls, MPI_BYTE, comm);
        st.rounds = 1;
        st.direct = 1;
    }
    else
    {
        int64_t cap = xchg_chunk_cap(n, unit, budget);
        int64_t local_rounds = 0, rounds, sstage = 0, rstage = 0;
        for (int p = 0; p < n; p++)
        {
            int64_t a = xchg_rounds_for(send_bytes[p], cap);
            int64_t b = xchg_rounds_for(recv_bytes[p], cap);
            if (a > local_rounds) local_rounds = a;
            if (b > local_rounds) local_rounds = b;
            sstage += (send_bytes[p] < cap) ? send_bytes[p] : cap;
            rstage += (recv_bytes[p] < cap) ? recv_bytes[p] : cap;
        }
        MPI_Allreduce(&local_rounds, &rounds, 1, MPI_INT64_T, MPI_MAX, comm);

        uint8_t *sbuf = malloc(sstage > 0 ? sstage : 1);
        uint8_t *rbuf = malloc(rstage > 0 ? rstage : 1);
        if (sbuf == NULL || rbuf == NULL) die(comm, "out of memory for staging", sstage + rstage);

        for (int64_t k = 0; k < rounds; k++)
        {
            int64_t s = xchg_plan_round(send_bytes, n, unit, cap, k, scounts, sdispls);
            int64_t r = xchg_plan_round(recv_bytes, n, unit, cap, k, rcounts, rdispls);
            if (s < 0 || r < 0) die(comm, "round planning failed", s < 0 ? s : r);
            for (int p = 0; p < n; p++)
            {
                if (scounts[p] == 0) continue;
                int64_t off = k * cap;
                if (src->buf)
                    memcpy(sbuf + sdispls[p], src->buf + sprefix[p] + off, scounts[p]);
                else
                    src->fill(src->ctx, p, off, scounts[p], sbuf + sdispls[p]);
            }
            MPI_Alltoallv(sbuf, scounts, sdispls, MPI_BYTE,
                          rbuf, rcounts, rdispls, MPI_BYTE, comm);
            for (int p = 0; p < n; p++)
            {
                if (rcounts[p] == 0) continue;
                int64_t off = k * cap;
                if (dst->buf)
                    memcpy(dst->buf + rprefix[p] + off, rbuf + rdispls[p], rcounts[p]);
                else
                    dst->drain(dst->ctx, p, off, rcounts[p], rbuf + rdispls[p]);
            }
        }
        free(sbuf);
        free(rbuf);
        st.rounds = rounds;
        st.cap = cap;
        st.staging_bytes = sstage + rstage;
    }

    free(sprefix); free(rprefix);
    free(scounts); free(sdispls); free(rcounts); free(rdispls);
    if (stats) *stats = st;
}
