#!/bin/bash
# Phase 2 regression check: build pga-sort, run it on yeast7, sort the output
# lexicographically (Phase 1 convention), and confirm it is byte-identical to
# tests/fixtures/yeast/gold.tuples.
#
# Validates: extract -> FastGA-format record pack -> msd_sort -> record-to-text decode.
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

echo "[1/4] Building pga-sort..."
( cd "$ROOT" && make pga-sort >/dev/null )

echo "[2/4] Extracting + sorting records via msd_sort..."
"$ROOT/pga-sort" "$GDB_STEM" "$WORK/our.msd-sorted.tuples"

echo "[3/4] Lex-sorting (msd_sort orders by k-mer; gold is sorted lexicographically)..."
sort -S 4G --parallel=4 "$WORK/our.msd-sorted.tuples" > "$WORK/our.tuples"

echo "[4/4] Comparing against gold.tuples..."
if cmp -s "$WORK/our.tuples" "$GOLD"; then
    echo "PASS: pga-sort output matches $GOLD after lexicographic sort"
    exit 0
else
    echo "FAIL: pga-sort output differs from $GOLD" >&2
    echo "  ours:  $(wc -l < "$WORK/our.tuples") lines" >&2
    echo "  gold:  $(wc -l < "$GOLD") lines" >&2
    diff <(head -10 "$WORK/our.tuples") <(head -10 "$GOLD") | head -20 >&2
    exit 1
fi
