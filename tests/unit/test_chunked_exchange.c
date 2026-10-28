// ============================================================================
// R-06 -- Exchanges larger than INT_MAX bytes.  See tests/CATALOG.md.
// ============================================================================
//
// THE PROBLEM THIS GUARDS.  Both MPI exchanges (stage-2 k-mer records in
// run_stage2, seeds in redistribute_seeds) used one MPI_Alltoallv with int
// counts and displacements and aborted once a rank sent or received more than
// INT_MAX bytes. primates4.chr6 needs ~8 GB of records and ~23 GB of seeds, so
// it could not run at small rank counts. chunked_alltoallv (src/chunked_exchange.c)
// moves the data in rounds of at most `cap` bytes per peer.
//
// The planner is MPI-free, so this tests multi-gigabyte sizes without
// allocating them, and then simulates a whole multi-rank exchange in-process to
// show the chunked result is byte-identical to one unbounded Alltoallv.
// The MPI path itself is exercised by the e2e/stage-2 runs with a tiny
// PGA_EXCHANGE_ROUND_BYTES, which forces many rounds.

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "chunked_exchange.h"

#define GiB ((int64_t) 1 << 30)

static void r06_cap_bounds(void)
{
    int64_t cap;

    cap = xchg_chunk_cap(8, 16, GiB);
    CHECK_EQ_I(cap, GiB / 8);
    CHECK(cap % 16 == 0);

    cap = xchg_chunk_cap(1, 16, (int64_t) 1 << 40);    // budget far above INT_MAX
    CHECK(cap <= INT_MAX && cap % 16 == 0 && cap > INT_MAX - 16);

    cap = xchg_chunk_cap(1024, 9, 4096);                // budget smaller than npeers*unit
    CHECK_EQ_I(cap, 9);

    cap = xchg_chunk_cap(3, 9, 1000);                   // 333 = 37 seeds
    CHECK_EQ_I(cap, 333);

    for (int n = 1; n <= 4096; n *= 2) {
        cap = xchg_chunk_cap(n, 16, 4 * GiB);
        CHECK_MSG(cap * n <= INT_MAX, "n=%d cap=%lld: round could exceed INT_MAX", n, (long long) cap);
    }
}

// Walk every round of a plan and check it is valid and complete.
static void check_plan(const int64_t *bytes, int n, int unit, int64_t cap, const char *label)
{
    int64_t rounds = 0, *moved = calloc(n, sizeof(int64_t));
    int *counts = malloc(n * sizeof(int)), *displs = malloc(n * sizeof(int));

    for (int p = 0; p < n; p++)
        if (xchg_rounds_for(bytes[p], cap) > rounds) rounds = xchg_rounds_for(bytes[p], cap);

    for (int64_t k = 0; k < rounds; k++) {
        int64_t tot = xchg_plan_round(bytes, n, unit, cap, k, counts, displs);
        CHECK_MSG(tot >= 0 && tot <= INT_MAX, "%s round %lld: total %lld", label,
                  (long long) k, (long long) tot);
        int64_t run = 0;
        for (int p = 0; p < n; p++) {
            CHECK_MSG(counts[p] % unit == 0, "%s: split record (peer %d round %lld)",
                      label, p, (long long) k);
            CHECK(counts[p] <= cap);
            CHECK(displs[p] == run);
            // chunks of one peer are contiguous: this one starts where the last ended
            if (counts[p] > 0) CHECK(moved[p] == k * cap);
            moved[p] += counts[p];
            run += counts[p];
        }
    }
    for (int p = 0; p < n; p++)
        CHECK_MSG(moved[p] == bytes[p], "%s peer %d: moved %lld of %lld", label, p,
                  (long long) moved[p], (long long) bytes[p]);
    // one round past the end moves nothing
    CHECK(xchg_plan_round(bytes, n, unit, cap, rounds, counts, displs) == 0);
    free(moved); free(counts); free(displs);
}

static void r06_plans_beyond_int_max(void)
{
    int64_t a[4] = { 5 * GiB, 0, 16, 3 * GiB + 32 };
    check_plan(a, 4, 16, xchg_chunk_cap(4, 16, GiB), "mixed");

    int64_t skew[4] = { 10 * GiB, 16, 16, 16 };       // one peer dominates
    int64_t cap = xchg_chunk_cap(4, 16, GiB);
    check_plan(skew, 4, 16, cap, "skew");
    CHECK_EQ_I(xchg_rounds_for(skew[0], cap), 40);
    CHECK_EQ_I(xchg_round_count(skew[1], cap, 1), 0);

    int64_t seeds[3] = { 23 * GiB / 9 * 9, 9, 7 * GiB / 9 * 9 };   // 9-byte seeds
    check_plan(seeds, 3, 9, xchg_chunk_cap(3, 9, GiB), "seeds");

    int64_t one[1] = { 4 * GiB };                       // a single rank to itself
    check_plan(one, 1, 16, xchg_chunk_cap(1, 16, GiB), "self");
}

static void r06_rejects_bad_sizes(void)
{
    int counts[2], displs[2];
    int64_t neg[2] = { -16, 16 }, odd[2] = { 16, 17 };
    CHECK_EQ_I(xchg_plan_round(neg, 2, 16, 64, 0, counts, displs), -1);
    CHECK_EQ_I(xchg_plan_round(odd, 2, 16, 64, 0, counts, displs), -2);
}

// In-process model of an N-rank exchange: rank s sends size[s][r] bytes to r.
// Reference = one unbounded Alltoallv (recv[r] = concat over s of s's block for r).
// Chunked = the same round plan chunked_alltoallv uses, placing each chunk at
// rprefix + k*cap. Must be byte-identical.
static void simulate(int n, int unit, int64_t budget, unsigned seed, const char *label)
{
    int64_t size[8][8], sprefix[8][8], rprefix[8][8], stot[8] = { 0 }, rtot[8] = { 0 };
    uint8_t *send[8], *ref[8], *got[8];
    srand(seed);
    for (int s = 0; s < n; s++)
        for (int r = 0; r < n; r++) {
            int64_t recs = (rand() % 5 == 0) ? 0 : rand() % 200;   // some empty, uneven
            if (s == 0 && r == n - 1) recs = 1500;                  // one heavy pair
            size[s][r] = recs * unit;
        }
    for (int s = 0; s < n; s++) {
        for (int r = 0; r < n; r++) { sprefix[s][r] = stot[s]; stot[s] += size[s][r]; }
        send[s] = malloc(stot[s] + 1);
        for (int64_t i = 0; i < stot[s]; i++) send[s][i] = (uint8_t) rand();
    }
    for (int r = 0; r < n; r++) {
        for (int s = 0; s < n; s++) { rprefix[r][s] = rtot[r]; rtot[r] += size[s][r]; }
        ref[r] = malloc(rtot[r] + 1);
        got[r] = calloc(rtot[r] + 1, 1);
        for (int s = 0; s < n; s++)
            memcpy(ref[r] + rprefix[r][s], send[s] + sprefix[s][r], size[s][r]);
    }

    int64_t cap = xchg_chunk_cap(n, unit, budget), rounds = 0;
    for (int s = 0; s < n; s++)
        for (int r = 0; r < n; r++)
            if (xchg_rounds_for(size[s][r], cap) > rounds) rounds = xchg_rounds_for(size[s][r], cap);
    for (int64_t k = 0; k < rounds; k++)
        for (int s = 0; s < n; s++)
            for (int r = 0; r < n; r++) {
                int64_t c = xchg_round_count(size[s][r], cap, k);
                if (c) memcpy(got[r] + rprefix[r][s] + k * cap, send[s] + sprefix[s][r] + k * cap, c);
            }
    for (int r = 0; r < n; r++) {
        CHECK_MSG(memcmp(ref[r], got[r], rtot[r]) == 0,
                  "%s: rank %d receive buffer differs from single Alltoallv", label, r);
        free(ref[r]); free(got[r]);
    }
    for (int s = 0; s < n; s++) free(send[s]);
    CHECK_MSG(rounds > 1, "%s: expected several rounds, got %lld", label, (long long) rounds);
}

static void r06_simulated_exchange_is_exact(void)
{
    simulate(5, 16, 5 * 16 * 3, 1u, "records, cap 3 records");
    simulate(8, 9, 8 * 9 * 7, 2u, "seeds, cap 7 seeds");
    simulate(3, 16, 3 * 16, 3u, "cap 1 record");
    simulate(1, 9, 9 * 4, 4u, "single rank");
}

static void r06_budget_env(void)
{
    unsetenv("PGA_EXCHANGE_ROUND_BYTES");
    CHECK_EQ_I(xchg_round_budget_from_env(), XCHG_DEFAULT_ROUND_BYTES);
    setenv("PGA_EXCHANGE_ROUND_BYTES", "4096", 1);
    CHECK_EQ_I(xchg_round_budget_from_env(), 4096);
    setenv("PGA_EXCHANGE_ROUND_BYTES", "0", 1);
    CHECK_EQ_I(xchg_round_budget_from_env(), XCHG_DEFAULT_ROUND_BYTES);
    setenv("PGA_EXCHANGE_ROUND_BYTES", "junk", 1);
    CHECK_EQ_I(xchg_round_budget_from_env(), XCHG_DEFAULT_ROUND_BYTES);
    unsetenv("PGA_EXCHANGE_ROUND_BYTES");
}

int main(void)
{
    printf("test_chunked_exchange\n");
    RUN(r06_cap_bounds);
    RUN(r06_plans_beyond_int_max);
    RUN(r06_rejects_bad_sizes);
    RUN(r06_simulated_exchange_is_exact);
    RUN(r06_budget_env);
    return t_report("test_chunked_exchange");
}
