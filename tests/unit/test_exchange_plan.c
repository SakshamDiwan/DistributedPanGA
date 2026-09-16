// ============================================================================
// R-02 -- The MPI_Alltoallv count/displacement planner.  See tests/CATALOG.md.
// ============================================================================
//
// THE PROBLEM THIS GUARDS.  MPI_Alltoallv takes its counts and displacements
// as `int`.  A displacement is a running offset into the send buffer, so it
// grows with the TOTAL bytes being sent, not with any one destination's share.
// That makes two separate overflow risks:
//
//   per-slot     one destination is handed more than INT_MAX bytes;
//   cumulative   every destination is individually fine, but their SUM passes
//                INT_MAX, so a later displacement wraps negative.
//
// The original code in redistribute_seeds checked the first and missed the
// second, so three destinations of 1 GiB each produced
// displs = {0, 1073741824, -2147483648} and MPI would have read outside the
// send buffer.  PR-2 moved this arithmetic into src/seed_exchange.c so it
// could be tested, and the cumulative guard was added there.
//
// HOW THIS TESTS GIGABYTE CASES IN MICROSECONDS: plan_byte_exchange takes the
// sizes as an array of integers and never touches a buffer, so the test
// describes a 3 GiB exchange without allocating a single byte.  That is also
// why the helper is MPI-free and returns a code instead of calling MPI_Abort:
// the caller keeps the abort, and the arithmetic stays testable in isolation.

#include <stdio.h>
#include <limits.h>
#include <stdint.h>

#include "check.h"
#include "seed_exchange.h"

// The ordinary case: counts pass through unchanged and displacements are their
// running prefix sums, so destination r's bytes begin at displs[r].
static void r02_ok_small(void)
{
    int64_t sizes[3] = { 100, 200, 300 };
    int counts[3], displs[3];

    CHECK_EQ_I(plan_byte_exchange(sizes, 3, counts, displs), 0);
    CHECK_EQ_I(counts[0], 100); CHECK_EQ_I(counts[1], 200); CHECK_EQ_I(counts[2], 300);
    CHECK_EQ_I(displs[0], 0);   CHECK_EQ_I(displs[1], 100); CHECK_EQ_I(displs[2], 300);
}

// Boundary: a cumulative total of exactly INT_MAX is still representable and
// must be ACCEPTED.  Tested because an off-by-one in the guard would reject
// valid work, which is a silent capacity regression rather than a crash.
static void r02_ok_at_limit(void)
{
    int64_t sizes[2] = { (int64_t) INT_MAX - 10, 10 };
    int counts[2], displs[2];

    CHECK_EQ_I(plan_byte_exchange(sizes, 2, counts, displs), 0);
    CHECK_EQ_I(displs[1], INT_MAX - 10);
}

// The per-slot guard, which existed before PR-2: one destination alone exceeds
// INT_MAX.  Kept so the refactor is pinned as behaviour-preserving.
static void r02_reject_single_slot(void)
{
    int64_t sizes[2] = { (int64_t) INT_MAX + 1, 5 };
    int counts[2], displs[2];

    CHECK_MSG(plan_byte_exchange(sizes, 2, counts, displs) != 0,
              "a single destination above INT_MAX must be rejected");
}

static void r02_reject_negative(void)
{
    int64_t sizes[2] = { 10, -1 };
    int counts[2], displs[2];
    CHECK(plan_byte_exchange(sizes, 2, counts, displs) != 0);
}

// The case the original code missed: every individual slot fits comfortably,
// and only the RUNNING TOTAL overflows.  3 x 1 GiB = 3 GiB, so displs[2] would
// be (int) 2147483648, i.e. negative.  The first two assertions state the
// setup explicitly -- each slot under INT_MAX, the sum over it -- so the
// fixture cannot drift into testing the per-slot guard by accident.
//
// The second half uses eight mid-sized destinations instead of three large
// ones, which is closer to what a real many-rank run looks like: no single
// destination is remotely near the limit, yet the total still is.
static void r02_reject_cumulative(void)
{
    int64_t sizes[3] = { 1 << 30, 1 << 30, 1 << 30 };
    int counts[3], displs[3];
    int rc;

    // each slot is individually fine
    CHECK(sizes[0] < (int64_t) INT_MAX);
    CHECK(sizes[0] * 3 > (int64_t) INT_MAX);

    rc = plan_byte_exchange(sizes, 3, counts, displs);
    CHECK_MSG(rc != 0, "3 x 1 GiB: cumulative total exceeds INT_MAX and must be "
                       "rejected (displs[2] would wrap negative)");

    // A wider case: many mid-sized destinations.
    {
        int64_t many[8];
        int c8[8], d8[8];
        int i, rc8;
        for (i = 0; i < 8; i++) many[i] = 400 * 1024 * 1024;   /* 400 MiB each */
        rc8 = plan_byte_exchange(many, 8, c8, d8);
        CHECK_MSG(rc8 != 0, "8 x 400 MiB: cumulative total exceeds INT_MAX");
    }
}

int main(void)
{
    printf("test_exchange_plan\n");
    RUN(r02_ok_small);
    RUN(r02_ok_at_limit);
    RUN(r02_reject_single_slot);
    RUN(r02_reject_negative);
    RUN(r02_reject_cumulative);
    return t_report("test_exchange_plan");
}
