#!/bin/bash
# Generate the small synthetic end-to-end fixture and its INDEPENDENT gold
# reference, produced by upstream FastGA rather than by this pipeline.
#
# Everything it writes lands in tests/fixtures/generated/, which .gitignore
# excludes: only the generator and this script are tracked, so a fresh clone
# reproduces the fixture rather than carrying binary data.
#
# Provenance (seed, tool revisions, parameters, hashes) is recorded in
# tests/fixtures/generated/MANIFEST.txt. Any result quoted from this gold set
# is meaningless without it.
#
# Usage:  scripts/make_gold_small.sh [seed]
# Env:    FASTGA_BIN  directory holding FAtoGDB/GIXmake/GIXshow/FastGA/ALNtoPAF
#                     (default $HOME/FASTGA)
#         FASTGA_SRC  optional FastGA source checkout, for recording its revision

set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=$ROOT/tests/fixtures/generated
GEN=$ROOT/tests/fixtures/small/gen_fixture.py
SEED=${1:-20260912}
STEM=small6
FASTGA_BIN=${FASTGA_BIN:-$HOME/FASTGA}

# Unique scratch, removed on exit. We only ever delete paths we created.
WORK=$(mktemp -d "${TMPDIR:-/tmp}/pga-gold-XXXXXX")
cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

# ---- prerequisites --------------------------------------------------------
missing=0
for t in FAtoGDB GIXmake GIXshow FastGA ALNtoPAF; do
    if [ ! -x "$FASTGA_BIN/$t" ]; then
        echo "ERROR: $FASTGA_BIN/$t not found or not executable" >&2
        missing=1
    fi
done
[ -x "$ROOT/notes/GIXmake_dump" ] || { echo "ERROR: notes/GIXmake_dump missing" >&2; missing=1; }
[ -f "$GEN" ] || { echo "ERROR: $GEN missing" >&2; missing=1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not found" >&2; missing=1; }
command -v gawk   >/dev/null || { echo "ERROR: gawk not found" >&2; missing=1; }
[ "$missing" -eq 0 ] || { echo "Install FastGA and set FASTGA_BIN, then retry." >&2; exit 1; }

export PATH="$FASTGA_BIN:$PATH"
mkdir -p "$OUT"

echo "[1/6] generating fixture (seed $SEED)"
python3 "$GEN" "$OUT/$STEM.fa" "$SEED"

echo "[2/6] FAtoGDB + GIXmake"
cd "$OUT"
rm -f "$STEM.1gdb"; rm -rf "$STEM.gix"
FAtoGDB "$STEM.fa"
GIXmake -T4 "$STEM" >/dev/null

echo "[3/6] gold.paf via FastGA self-alignment"
# ONE genome argument: self-alignment mode. "$STEM $STEM" would be cross-mode
# and produce a different alignment set.
rm -f gold.1aln gold.paf
FastGA -T4 -1:gold "$STEM" >/dev/null 2>&1
ALNtoPAF gold.1aln > gold.paf

echo "[4/6] gold.tuples via the instrumented GIXmake (rebuilds .gix as a side effect)"
rm -rf "$STEM.gix"
GIXMAKE_DUMP_PATH="$WORK/dump.tuples" "$ROOT/notes/GIXmake_dump" -T4 "$STEM" >/dev/null 2>&1
sort -S 1G "$WORK/dump.tuples" > gold.tuples

echo "[5/6] cross-check gold.tuples against GIXshow"
gcc -O2 -o "$WORK/normalize_gixshow" "$ROOT/scripts/normalize_gixshow.c"
GIXshow "$STEM.gix" | "$WORK/normalize_gixshow" | sort -S 1G > "$WORK/fastga_full.tuples"
if cmp -s gold.tuples "$WORK/fastga_full.tuples"; then
    echo "      gold.tuples is byte-identical to FastGA's .gix contents"
else
    echo "ERROR: gold.tuples differs from FastGA's .gix contents" >&2
    exit 1
fi

echo "[6/6] simple-LCP oracle"
GIXSHOW="$FASTGA_BIN/GIXshow" bash "$ROOT/scripts/build_lcp_oracle_simple.sh" \
    "$STEM.gix" gold.lcp.kmers.simple

# ---- provenance -----------------------------------------------------------
fastga_rev="unknown"
if [ -n "${FASTGA_SRC:-}" ] && [ -d "${FASTGA_SRC}/.git" ]; then
    fastga_rev=$(git -C "$FASTGA_SRC" rev-parse HEAD 2>/dev/null || echo unknown)
fi
{
    echo "# Provenance for the small end-to-end fixture and its gold reference."
    echo "# Regenerate with: scripts/make_gold_small.sh $SEED"
    echo
    echo "generated_utc:   $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "generator:       tests/fixtures/small/gen_fixture.py"
    echo "generator_seed:  $SEED"
    echo "fixture_params:  6 contigs x 8000 bp, two near-repeat families, 2% divergence"
    echo "kmer_params:     KMER=40 TMER=12 SMER=8 (FastGA defaults)"
    echo "fastga_bin:      $FASTGA_BIN"
    echo "fastga_revision: $fastga_rev"
    echo "pipeline_commit: $(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "pipeline_dirty:  $(test -n "$(git -C "$ROOT" status --porcelain)" && echo yes || echo no)"
    echo
    echo "alignments_in_gold_paf: $(wc -l < gold.paf)"
    echo
    echo "# sha256 of fixture, gold and the FastGA binaries that produced them"
    sha256sum "$STEM.fa" "$STEM.1gdb" gold.paf gold.tuples gold.lcp.kmers.simple 2>/dev/null
    for t in FAtoGDB GIXmake GIXshow FastGA ALNtoPAF; do sha256sum "$FASTGA_BIN/$t"; done
} > MANIFEST.txt

echo
echo "Done. Artifacts in $OUT (gitignored):"
ls -1 "$STEM".* gold.* MANIFEST.txt | sed 's/^/  /'
echo
echo "Gold contains $(wc -l < gold.paf) alignments. That number describes the"
echo "FIXTURE, not this pipeline -- validation is the comparison in"
echo "scripts/verify_e2e_small.sh, not the count."
