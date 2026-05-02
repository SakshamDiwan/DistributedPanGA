#!/bin/bash
set -euo pipefail

ROOT=$HOME/DistributedPanGA
FIXTURES=$ROOT/tests/fixtures/yeast
INPUT=$FIXTURES/yeast7

if [ ! -f "$INPUT.fa.gz" ]; then
    echo "ERROR: $INPUT.fa.gz not found." >&2
    exit 1
fi

cd "$FIXTURES"

echo "[1/5] Building GDB and GIX from yeast7.fa.gz..."
rm -f yeast7.1gdb yeast7.gix
~/FASTGA/FAtoGDB yeast7.fa.gz
~/FASTGA/GIXmake yeast7

echo "[2/5] gold.kmers via GIXshow..."
~/FASTGA/GIXshow yeast7.gix | sort -S 4G --parallel=4 > gold.kmers

echo "[3/5] gold.paf via FastGA self-alignment..."
~/FASTGA/FastGA -1:gold yeast7 yeast7
~/FASTGA/ALNtoPAF gold.1aln > gold.paf

echo "[4/5] gold.tuples via instrumented GIXmake..."
# The instrumented run rebuilds yeast7.gix as a side effect — that's fine,
# the new .gix is identical to the one from step 1.
rm -f yeast7.gix
GIXMAKE_DUMP_PATH=/tmp/yeast_dump.tuples \
    "$ROOT/notes/GIXmake_dump" -v yeast7
sort -S 4G --parallel=4 /tmp/yeast_dump.tuples > gold.tuples
rm /tmp/yeast_dump.tuples

echo "[5/5] Verifying byte-identical match against FastGA..."
if [ ! -f /tmp/normalize_gixshow ]; then
    gcc -O2 -o /tmp/normalize_gixshow "$ROOT/scripts/normalize_gixshow.c"
fi
~/FASTGA/GIXshow yeast7.gix | /tmp/normalize_gixshow \
    | sort -S 4G --parallel=4 > /tmp/fastga_full.tuples
if cmp -s gold.tuples /tmp/fastga_full.tuples; then
    echo "  gold.tuples is byte-identical to FastGA's .gix contents ✓"
else
    echo "  WARNING: gold.tuples DIFFERS from FastGA's .gix contents" >&2
    diff <(head -5 gold.tuples) <(head -5 /tmp/fastga_full.tuples) >&2
    exit 1
fi
rm -f /tmp/fastga_full.tuples

echo
echo "Done. Artifacts in $FIXTURES:"
ls -la gold.* yeast7.*
