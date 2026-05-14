#!/bin/bash
# Phase 3 stage 2 regression: build pga-mpi, run at N ∈ {1,2,4,8} ranks,
# concatenate per-rank outputs, and run two checks:
#   Tier 1 — multiset equality:  sort combined | cmp - gold
#   Tier 2 — global k-mer order on the *unsorted* concatenation
#            (proves the distributed sort actually distributed-sorted).
# Both checks must pass at every rank count.
#
# REQUIRES: an active Slurm allocation (salloc / sbatch). Uses srun.
# Exits 0 on success, non-zero on any failure.
#
# See implementation_plan.md §7a.7 for why both tiers are mandatory.

set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
FIXTURES=$ROOT/tests/fixtures/yeast
GOLD=$FIXTURES/gold.tuples
GDB_STEM=$FIXTURES/yeast7
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

if [ ! -f "$GOLD" ]; then
    echo "ERROR: $GOLD not found. Run scripts/make_gold.sh first." >&2
    exit 1
fi
if [ ! -f "$GDB_STEM.1gdb" ]; then
    echo "ERROR: $GDB_STEM.1gdb not found." >&2
    exit 1
fi
if ! command -v srun >/dev/null 2>&1; then
    echo "ERROR: srun not found. Run inside a Slurm allocation." >&2
    exit 1
fi

echo "[build] pga-mpi..."
( cd "$ROOT" && make pga-mpi >/dev/null )

RANK_COUNTS="${RANK_COUNTS:-1 2 4 8}"
overall_status=0

# Tier 2 walker: read concatenated tuple file, assert each line's k-mer (4th
# whitespace field, 20 hex chars) is >= the previous line's k-mer. Awk wins
# on memory because it streams.
tier2_check() {
    local file="$1"
    awk '
        NR == 1 { prev = $4; next }
        {
            if ($4 < prev) {
                printf("Tier 2 FAIL at line %d: kmer %s < prev %s\n", NR, $4, prev) > "/dev/stderr"
                exit 1
            }
            prev = $4
        }
    ' "$file"
}

for N in $RANK_COUNTS; do
    echo "=== N=$N ==="
    OUT_STEM="$WORK/out.N$N"
    COMBINED="$WORK/combined.N$N.tuples"
    SORTED="$WORK/sorted.N$N.tuples"

    echo "  [run] srun -n $N pga-mpi ..."
    srun -n "$N" "$ROOT/pga-mpi" "$GDB_STEM" "$OUT_STEM"

    echo "  [cat] concatenating ${N} per-rank outputs in rank order..."
    : > "$COMBINED"
    for ((r = 0; r < N; r++)); do
        cat "$OUT_STEM.$r.tuples" >> "$COMBINED"
    done

    echo "  [tier2] global k-mer order on unsorted concatenation..."
    if tier2_check "$COMBINED"; then
        echo "  Tier 2 PASS"
    else
        echo "  Tier 2 FAIL: cross-rank order is wrong" >&2
        overall_status=1
        continue
    fi

    echo "  [tier1] multiset equality vs gold (sort + cmp)..."
    sort -S 4G --parallel=4 "$COMBINED" > "$SORTED"
    if cmp -s "$SORTED" "$GOLD"; then
        echo "  Tier 1 PASS: N=$N matches gold"
    else
        echo "  Tier 1 FAIL: N=$N differs from gold" >&2
        echo "    ours: $(wc -l < "$SORTED") lines" >&2
        echo "    gold: $(wc -l < "$GOLD") lines" >&2
        diff <(head -10 "$SORTED") <(head -10 "$GOLD") | head -20 >&2 || true
        overall_status=1
    fi
done

if [ "$overall_status" -eq 0 ]; then
    echo ""
    echo "ALL PASS: stage 2 verified at N ∈ { $RANK_COUNTS } (Tier 1 + Tier 2)"
else
    echo ""
    echo "SOME FAILED" >&2
fi
exit "$overall_status"
