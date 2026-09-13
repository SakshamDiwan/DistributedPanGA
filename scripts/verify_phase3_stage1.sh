#!/bin/bash
# Phase 3 stage 1 regression: build pga-mpi-stage1, run it at N ∈ {1,2,4,8} ranks,
# and confirm the rank-0 sorted output matches gold.tuples (after external sort,
# same convention as verify_phase2.sh).
#
# REQUIRES: an active Slurm allocation (salloc / sbatch). Uses srun.
# Exits 0 on success, non-zero on any failure.

set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
# Overridable so the suite can run against a small synthetic fixture as well
# as the yeast reference. Defaults are the historical values.
FIXTURES=${PGA_FIXTURES:-$ROOT/tests/fixtures/yeast}
GOLD=$FIXTURES/gold.tuples
GDB_STEM=${PGA_GDB_STEM:-$FIXTURES/yeast7}
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

echo "[build] pga-mpi-stage1..."
( cd "$ROOT" && make pga-mpi-stage1 >/dev/null )

RANK_COUNTS="${RANK_COUNTS:-1 2 4 8}"
overall_status=0

for N in $RANK_COUNTS; do
    echo "=== N=$N ==="
    OUT="$WORK/out.N$N.tuples"
    SORTED="$WORK/sorted.N$N.tuples"

    echo "  [run] srun -n $N pga-mpi-stage1 ..."
    srun -n "$N" "$ROOT/pga-mpi-stage1" "$GDB_STEM" "$OUT"

    echo "  [sort] lex-sorting output..."
    sort -S 4G --parallel=4 "$OUT" > "$SORTED"

    echo "  [diff] cmp against gold.tuples..."
    if cmp -s "$SORTED" "$GOLD"; then
        echo "  PASS: N=$N matches gold"
    else
        echo "  FAIL: N=$N differs from gold" >&2
        echo "    ours: $(wc -l < "$SORTED") lines" >&2
        echo "    gold: $(wc -l < "$GOLD") lines" >&2
        diff <(head -10 "$SORTED") <(head -10 "$GOLD") | head -20 >&2 || true
        overall_status=1
    fi
done

if [ "$overall_status" -eq 0 ]; then
    echo ""
    echo "ALL PASS: stage 1 verified at N ∈ { $RANK_COUNTS }"
else
    echo ""
    echo "SOME FAILED" >&2
fi
exit "$overall_status"
