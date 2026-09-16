// ============================================================================
// P3 -- Deterministic contig-pair assignment.  Tests A-01..A-05.
// See tests/CATALOG.md.
// ============================================================================
//
// WHAT THIS DECIDES.  In self-alignment every contig is compared against every
// other, so for N contigs there are N*(N+1)/2 unordered pairs (i,j) with i<=j.
// Each pair must be aligned by exactly one rank.  contig_assignment_build
// hands out those pairs.
//
// THE POLICY (documented in contig_assignment.h):
//   1. weight each pair by clen[i] * clen[j], a proxy for how much work it is;
//   2. sort pairs by weight DESCENDING, breaking ties by (i,j) ASCENDING;
//   3. walk that order, giving each pair to the currently least-loaded rank,
//      breaking load ties by lower rank id.
//
// WHY DETERMINISM IS A CORRECTNESS PROPERTY, not a nicety: every rank runs
// this independently and they never compare answers.  The only reason that is
// safe is that identical inputs provably give identical output.  If two ranks
// disagreed about who owns a pair, it would either be aligned twice or not at
// all.  A-04 pins that directly, and proposal item O3 -- which replaces the
// length-product weights with measured seed densities -- must preserve it.
//
// The expected rank assignments in A-01 and A-02 were worked out by hand from
// the policy above, tracing the heap step by step.  They are not recorded
// output.
//
// Links only contig_assignment.c: it uses GDB struct fields but calls no GDB
// functions, so no lib/ sources are needed.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "contig_assignment.h"

static GDB *mk_gdb(const int64_t *clen, int n)
{
    GDB *g = calloc(1, sizeof(GDB));
    int  i;
    g->ncontig = n;
    g->contigs = calloc((size_t) n + 1, sizeof(GDB_CONTIG));
    for (i = 0; i < n; i++) {
        g->contigs[i].clen = clen[i];
        g->seqtot += clen[i];
        if (clen[i] > g->maxctg) g->maxctg = clen[i];
    }
    return g;
}

static void free_gdb(GDB *g) { free(g->contigs); free(g); }

// A-01 -- The smallest case where the tie-breaking rules are visible.
//
// Two contigs of equal length L on two ranks.  The three pairs (0,0), (0,1),
// (1,1) all weigh L*L, so the weight sort tells us nothing and every rule
// below it is exposed:
//
//   weights tie      -> order falls back to (i,j) ascending:
//                       (0,0), (0,1), (1,1)
//   (0,0) -> rank 0   (both ranks at load 0; lower rank id wins the tie)
//   (0,1) -> rank 1   (rank 0 now has L*L, rank 1 still 0)
//   (1,1) -> rank 0   (both at L*L again; lower rank id wins again)
//
// Expected: 0, 1, 0.
static void a01_two_equal(void)
{
    int64_t clen[2] = { 1000, 1000 };
    GDB *g = mk_gdb(clen, 2);
    ContigAssignment a;

    CHECK(contig_assignment_build(g, 2, &a) == 0);
    // all three weights tie at L^2 -> (i,j) ascending -> heap breaks by rank id
    CHECK_EQ_I(contig_assignment_owner(&a, 0, 0), 0);
    CHECK_EQ_I(contig_assignment_owner(&a, 0, 1), 1);
    CHECK_EQ_I(contig_assignment_owner(&a, 1, 1), 0);
    contig_assignment_free(&a);
    free_gdb(g);
}

// A-02 -- The same rules over six pairs, giving strict alternation.
//
// Three equal contigs make six pairs, all of equal weight, so again the order
// is purely (i,j) ascending:
//     (0,0) (0,1) (0,2) (1,1) (1,2) (2,2)
// and since every pair adds the same load, the least-loaded rank alternates
// each time:
//     0     1     0     1     0     1
//
// The alternation is a CONSEQUENCE of every weight tying, not a property of
// the algorithm -- with unequal contigs the pattern would not be regular.
// A-03 covers that case.
static void a02_three_equal(void)
{
    int64_t clen[3] = { 500, 500, 500 };
    GDB *g = mk_gdb(clen, 3);
    ContigAssignment a;
    int want[6] = { 0, 1, 0, 1, 0, 1 };
    int pair[6][2] = { {0,0},{0,1},{0,2},{1,1},{1,2},{2,2} };
    int i;

    CHECK(contig_assignment_build(g, 2, &a) == 0);
    for (i = 0; i < 6; i++)
        CHECK_MSG(contig_assignment_owner(&a, pair[i][0], pair[i][1]) == want[i],
                  "pair (%d,%d) -> %d, want %d", pair[i][0], pair[i][1],
                  contig_assignment_owner(&a, pair[i][0], pair[i][1]), want[i]);
    contig_assignment_free(&a);
    free_gdb(g);
}

// A-03 -- With irregular input, the structural guarantees still hold.
//
// Six contigs of distinct lengths over three ranks, where the exact placement
// is not worth hand-computing.  What must hold regardless:
//
//   COVERAGE    every pair i<=j has a valid owner in [0, world_size).  A pair
//               left unassigned would simply never be aligned.
//   TRIANGLE    entries with i>j are never written and stay -1.  The map is a
//               full N*N array but only the upper triangle is meaningful.
//   SYMMETRY    owner(i,j) == owner(j,i).  Callers index either way round, so
//               the accessor must normalise.
//   LOAD        the per-rank loads sum to the total of all pair weights, i.e.
//               work is neither invented nor dropped in the bookkeeping.
static void a03_coverage_and_load(void)
{
    int64_t clen[6] = { 900, 100, 550, 1200, 75, 640 };
    GDB *g = mk_gdb(clen, 6);
    ContigAssignment a;
    int64_t want_total = 0, got_total = 0;
    int i, j, r;

    CHECK(contig_assignment_build(g, 3, &a) == 0);

    for (i = 0; i < 6; i++)
        for (j = i; j < 6; j++) {
            int o = contig_assignment_owner(&a, i, j);
            CHECK_MSG(o >= 0 && o < 3, "pair (%d,%d) owner %d out of range", i, j, o);
            CHECK_MSG(contig_assignment_owner(&a, j, i) == o,
                      "owner not symmetric for (%d,%d)", i, j);
            want_total += clen[i] * clen[j];
        }

    // entries strictly below the diagonal are never written
    for (i = 0; i < 6; i++)
        for (j = 0; j < i; j++)
            CHECK_MSG(a.rank_for_pair[(int64_t) i * 6 + j] == -1,
                      "lower-triangle entry (%d,%d) = %d, want -1",
                      i, j, a.rank_for_pair[(int64_t) i * 6 + j]);

    for (r = 0; r < 3; r++) got_total += a.load_per_rank[r];
    CHECK_MSG(got_total == want_total, "load sum %lld, want %lld",
              (long long) got_total, (long long) want_total);

    contig_assignment_free(&a);
    free_gdb(g);
}

// A-04 -- Two builds from identical input are byte-identical.
//
// The whole design depends on this: ranks compute the map redundantly and
// never exchange it, so any nondeterminism (an unstable sort, a hash of a
// pointer, uninitialised padding) would cause ranks to silently disagree about
// pair ownership.  Note the deliberately duplicated contig lengths in the
// fixture -- ties are where an unstable sort would show itself.
static void a04_determinism(void)
{
    int64_t clen[7] = { 311, 97, 1024, 640, 640, 12, 888 };
    GDB *g = mk_gdb(clen, 7);
    ContigAssignment a, b;

    CHECK(contig_assignment_build(g, 4, &a) == 0);
    CHECK(contig_assignment_build(g, 4, &b) == 0);
    CHECK_MSG(memcmp(a.rank_for_pair, b.rank_for_pair, 49) == 0,
              "rank_for_pair differs between identical builds");
    CHECK_MSG(memcmp(a.load_per_rank, b.load_per_rank, 4 * sizeof(int64_t)) == 0,
              "load_per_rank differs between identical builds");
    contig_assignment_free(&a);
    contig_assignment_free(&b);
    free_gdb(g);
}

// A-05 -- The int8_t storage limit is enforced, and enforced in the right place.
//
// rank_for_pair is an int8_t array, so it cannot represent a rank id above
// 127.  The code rejects world_size > 127 rather than silently truncating,
// which would be far worse: pairs would be assigned to ranks that do not
// exist, or wrap to negative.
//
// Pinned from both sides -- 127 must still WORK, not just 128 must fail --
// because a guard that was off by one, or that rejected everything, would
// otherwise look correct.  Run under ASan so the rejection path is also
// checked for leaking the partially built map.
static void a05_world_size_bound(void)
{
    int64_t clen[4] = { 100, 200, 300, 400 };
    GDB *g = mk_gdb(clen, 4);
    ContigAssignment a;

    CHECK_MSG(contig_assignment_build(g, 127, &a) == 0, "world_size 127 should succeed");
    contig_assignment_free(&a);

    memset(&a, 0, sizeof a);
    CHECK_MSG(contig_assignment_build(g, 128, &a) != 0,
              "world_size 128 exceeds int8_t storage and must be rejected");
    free_gdb(g);
}

int main(void)
{
    printf("test_assign\n");
    RUN(a01_two_equal);
    RUN(a02_three_equal);
    RUN(a03_coverage_and_load);
    RUN(a04_determinism);
    RUN(a05_world_size_bound);
    return t_report("test_assign");
}
