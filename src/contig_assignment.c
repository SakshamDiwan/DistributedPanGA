#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "contig_assignment.h"

typedef struct {
    int     i;
    int     j;
    int64_t weight;
} Pair;

static int pair_cmp_weight_desc(const void *a, const void *b)
{
    const Pair *pa = (const Pair *) a;
    const Pair *pb = (const Pair *) b;
    // Sort by weight descending; tie-break by (i, j) ascending for
    // determinism across machines/libcs.
    if (pa->weight != pb->weight) return (pa->weight < pb->weight) ? 1 : -1;
    if (pa->i     != pb->i)      return pa->i - pb->i;
    return pa->j - pb->j;
}

// Min-heap over ranks keyed by (load, rank). Heap[0] is unused; entries
// start at heap[1]. Each entry is a rank id; load_per_rank holds the key.
typedef struct {
    int     *rank;
    int64_t *load;
    int      size;
} MinHeap;

static void heap_swap(MinHeap *h, int a, int b)
{
    int t = h->rank[a]; h->rank[a] = h->rank[b]; h->rank[b] = t;
}

// Heap key: (load[rank], rank). Lower load wins; ties broken by lower rank
// (so a later-assigned rank can't leapfrog an earlier one with same load).
static int heap_less(MinHeap *h, int a, int b)
{
    int ra = h->rank[a], rb = h->rank[b];
    if (h->load[ra] != h->load[rb]) return h->load[ra] < h->load[rb];
    return ra < rb;
}

static void heap_sift_down(MinHeap *h, int i)
{
    while (1) {
        int l = 2 * i, r = 2 * i + 1, smallest = i;
        if (l <= h->size && heap_less(h, l, smallest)) smallest = l;
        if (r <= h->size && heap_less(h, r, smallest)) smallest = r;
        if (smallest == i) return;
        heap_swap(h, i, smallest);
        i = smallest;
    }
}

int contig_assignment_build(const GDB *gdb, int world_size, ContigAssignment *out)
{
    if (world_size <= 0 || gdb == NULL || gdb->ncontig <= 0) return -1;

    int     N        = gdb->ncontig;
    int64_t npairs64 = (int64_t) N * (N + 1) / 2;
    if (npairs64 > (int64_t) 1 << 30) {
        // Sanity guard; with 32k contigs npairs is ~500M which is fine but
        // beyond that we'd want a different storage shape.
        fprintf(stderr, "contig_assignment_build: too many pairs (%lld)\n",
                (long long) npairs64);
        return -1;
    }
    int npairs = (int) npairs64;

    out->world_size    = world_size;
    out->ncontig       = N;
    out->rank_for_pair = (int8_t *)  malloc((size_t) N * N);
    out->load_per_rank = (int64_t *) calloc(world_size, sizeof(int64_t));
    if (out->rank_for_pair == NULL || out->load_per_rank == NULL) {
        free(out->rank_for_pair);
        free(out->load_per_rank);
        return -1;
    }
    if (world_size > 127) {
        fprintf(stderr, "contig_assignment_build: world_size %d > 127 "
                        "exceeds int8_t storage\n", world_size);
        free(out->rank_for_pair);
        free(out->load_per_rank);
        return -1;
    }
    memset(out->rank_for_pair, -1, (size_t) N * N);

    Pair *pairs = (Pair *) malloc((size_t) npairs * sizeof(Pair));
    if (pairs == NULL) {
        free(out->rank_for_pair);
        free(out->load_per_rank);
        return -1;
    }
    int p = 0;
    for (int i = 0; i < N; i++) {
        int64_t li = gdb->contigs[i].clen;
        for (int j = i; j < N; j++) {
            int64_t lj = gdb->contigs[j].clen;
            pairs[p].i      = i;
            pairs[p].j      = j;
            pairs[p].weight = li * lj;
            p++;
        }
    }
    qsort(pairs, npairs, sizeof(Pair), pair_cmp_weight_desc);

    MinHeap heap;
    heap.rank = (int *) malloc(sizeof(int) * (world_size + 1));
    heap.load = out->load_per_rank;
    heap.size = world_size;
    if (heap.rank == NULL) {
        free(pairs);
        free(out->rank_for_pair);
        free(out->load_per_rank);
        return -1;
    }
    for (int r = 0; r < world_size; r++) heap.rank[r + 1] = r;
    // Already a valid heap (all loads 0, ties broken by rank id ascending).

    for (int k = 0; k < npairs; k++) {
        int r = heap.rank[1];
        out->rank_for_pair[(int64_t) pairs[k].i * N + pairs[k].j] = (int8_t) r;
        out->load_per_rank[r] += pairs[k].weight;
        heap_sift_down(&heap, 1);
    }

    free(heap.rank);
    free(pairs);
    return 0;
}

void contig_assignment_free(ContigAssignment *a)
{
    if (a == NULL) return;
    free(a->rank_for_pair); a->rank_for_pair = NULL;
    free(a->load_per_rank); a->load_per_rank = NULL;
    a->ncontig    = 0;
    a->world_size = 0;
}

void contig_assignment_print(const ContigAssignment *a, const GDB *gdb)
{
    fprintf(stderr, "ContigAssignment: %d contigs, %d ranks\n",
            a->ncontig, a->world_size);
    int64_t min_load = a->load_per_rank[0];
    int64_t max_load = a->load_per_rank[0];
    int64_t sum_load = 0;
    int     n_pairs_per_rank[128] = {0};
    for (int r = 0; r < a->world_size; r++) {
        if (a->load_per_rank[r] < min_load) min_load = a->load_per_rank[r];
        if (a->load_per_rank[r] > max_load) max_load = a->load_per_rank[r];
        sum_load += a->load_per_rank[r];
    }
    for (int i = 0; i < a->ncontig; i++) {
        for (int j = i; j < a->ncontig; j++) {
            int8_t r = a->rank_for_pair[(int64_t) i * a->ncontig + j];
            if (r >= 0 && r < 128) n_pairs_per_rank[r]++;
        }
    }
    for (int r = 0; r < a->world_size; r++) {
        double pct = sum_load > 0
            ? 100.0 * (double) a->load_per_rank[r] / (double) sum_load
            : 0.0;
        fprintf(stderr, "  rank %2d: %d pairs, load %lld (%.1f%%)\n",
                r, n_pairs_per_rank[r], (long long) a->load_per_rank[r], pct);
    }
    double imbalance = (double) (max_load - min_load) / (double) (sum_load / a->world_size);
    fprintf(stderr, "  load imbalance: %.2f%% of mean\n", 100.0 * imbalance);
    (void) gdb;
}
