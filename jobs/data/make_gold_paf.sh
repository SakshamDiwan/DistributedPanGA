#!/bin/bash
# Build an independent FastGA self-alignment gold (PAF only) for one dataset.
# Usage: make_gold_paf.sh <dataset-dir> <stem> <source.fa.gz> <threads> <tag>
#   Creates <dir>/<stem>.fa.gz (symlink), <stem>.1gdb/.gix, <tag>.1aln, <tag>.paf
#   and appends provenance to <dir>/MANIFEST.txt. Phase-level gold (tuples/kmers)
#   is NOT built: it is ~80 B/base and only needed for yeast7's verifiers.
set -euo pipefail
DIR=$1 STEM=$2 SRC=$3 T=$4 TAG=$5
FG=$HOME/FASTGA
cd "$DIR"
[ -e "$STEM.fa.gz" ] || ln -s "$SRC" "$STEM.fa.gz"
mkdir -p tmp
tv() { /usr/bin/time -f "    %e s wall, %M kB maxRSS" "$@"; }
if [ ! -s "$STEM.gix" ]; then
  echo "  FAtoGDB";  rm -f "$STEM.1gdb"; tv $FG/FAtoGDB "$STEM.fa.gz"
  echo "  GIXmake -T$T"; tv $FG/GIXmake -T$T -Ptmp "$STEM"
fi
echo "  FastGA -1:$TAG -T$T"; rm -f "$TAG.1aln" "$TAG.paf"
tv $FG/FastGA -1:"$TAG" -T$T -Ptmp "$STEM"
echo "  ALNtoPAF"; tv $FG/ALNtoPAF -T$T "$TAG.1aln" > "$TAG.paf"
rmdir tmp 2>/dev/null || true
n=$(wc -l < "$TAG.paf")
echo "  $TAG.paf: $n alignments"
{
  echo "# $TAG  ($(date -u +%FT%TZ), slurm job ${SLURM_JOB_ID:-none}, host $(hostname))"
  echo "source:      $(readlink -f "$SRC")"
  echo "fastga_bin:  $FG (sha256 FastGA $(sha256sum $FG/FastGA | cut -c1-16)…)"
  echo "command:     FAtoGDB $STEM.fa.gz; GIXmake -T$T $STEM; FastGA -1:$TAG -T$T $STEM; ALNtoPAF $TAG.1aln > $TAG.paf"
  echo "defaults:    -f10 (FREQ), -k40, -c85 -s1000 -l100 -i.7"
  echo "alignments:  $n"
  sha256sum "$(readlink -f "$SRC")" "$STEM.1gdb" "$TAG.paf"
  echo
} >> MANIFEST.txt
