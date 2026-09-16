#!/bin/bash
# Phase 1 regression check: build pga-extract, run it on yeast7, sort the output,
# and confirm it is byte-identical to tests/fixtures/yeast/gold.tuples.
#
# Exits 0 on success, non-zero on any failure (build, run, missing fixtures, diff).

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

echo "[1/4] Building pga-extract..."
( cd "$ROOT" && make pga-extract >/dev/null )

echo "[2/4] Extracting tuples from yeast7..."
"$ROOT/pga-extract" "$GDB_STEM" "$WORK/our.unsorted.tuples"

echo "[3/4] Sorting..."
sort -S 4G --parallel=4 "$WORK/our.unsorted.tuples" > "$WORK/our.tuples"

echo "[4/4] Comparing against gold.tuples..."
if cmp -s "$WORK/our.tuples" "$GOLD"; then
    echo "PASS: pga-extract output is byte-identical to $GOLD"
    exit 0
else
    echo "FAIL: pga-extract output differs from $GOLD" >&2
    echo "  ours:  $(wc -l < "$WORK/our.tuples") lines" >&2
    echo "  gold:  $(wc -l < "$GOLD") lines" >&2
    echo "First divergence:" >&2
    diff <(head -200 "$WORK/our.tuples") <(head -200 "$GOLD") | head -20 >&2
    exit 1
fi
