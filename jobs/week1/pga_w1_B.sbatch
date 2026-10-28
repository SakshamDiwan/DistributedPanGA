#!/bin/bash
# Week 1, Job B: independent FastGA gold for yeast7 (two contig-order candidates).
#   yeast/        = cerevisiae.fa.gz as downloaded (S288C first)       -> full make_gold.sh
#   yeast-alpha/  = same 112 seqs, samples in the appendix's order     -> GDB + FastGA PAF only
# Whichever reproduces the historical 43,126 alignments identifies the input order.
# No pipeline binary is run here. Submit: sbatch ~/pga-jobs/week1/pga_w1_B.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --time=01:00:00
#SBATCH --job-name=pga-w1-B
#SBATCH --output=%x.%j.log
set -uo pipefail
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
R=$HOME/DistributedPanGA; F=$R/tests/fixtures
export PATH="$HOME/FASTGA:$PATH"
# Fixtures are symlinks into $SCRATCH/pga-fixtures (home quota is 40 GiB; gold.kmers alone is ~5.4 GB).
readlink -f $F/yeast $F/yeast-alpha
myquota 2>&1 | head -3
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)"
echo "PIPELINE COMMIT: $(git -C $R rev-parse HEAD)"
sha256sum $HOME/FASTGA/{FAtoGDB,GIXmake,GIXshow,FastGA,ALNtoPAF} $R/notes/GIXmake_dump
sha256sum $F/yeast/yeast7.fa.gz $F/yeast-alpha/yeast7.fa.gz

echo; echo "### 1. full gold, file-order candidate (scripts/make_gold.sh) ###"
t0=$(date +%s)
bash $R/scripts/make_gold.sh > $F/yeast/make_gold.log 2>&1; echo "make_gold rc=$?  ($(( $(date +%s)-t0 )) s)"
tail -15 $F/yeast/make_gold.log
( cd $F/yeast && wc -l gold.paf gold.tuples gold.lcp.kmers.simple 2>&1; sha256sum yeast7.1gdb gold.paf gold.tuples gold.lcp.kmers.simple 2>&1 )

echo; echo "### 2. FastGA PAF, alphabetical-order candidate ###"
cd $F/yeast-alpha && rm -f yeast7.1gdb yeast7.gix gold.1aln gold.paf
t0=$(date +%s)
FAtoGDB yeast7.fa.gz && GIXmake yeast7 && FastGA -1:gold yeast7 && ALNtoPAF gold.1aln > gold.paf
echo "alpha rc=$?  ($(( $(date +%s)-t0 )) s)"
wc -l gold.paf; sha256sum yeast7.1gdb gold.paf

echo; echo "### 3. order-independence check (PAF lines, sorted) ###"
for d in yeast yeast-alpha; do [ -s $F/$d/gold.paf ] && echo "$d: $(wc -l < $F/$d/gold.paf) alignments" || echo "$d: gold.paf MISSING/EMPTY"; done
if [ -s $F/yeast/gold.paf ] && [ -s $F/yeast-alpha/gold.paf ]; then
  cmp -s <(LC_ALL=C sort $F/yeast/gold.paf) <(LC_ALL=C sort $F/yeast-alpha/gold.paf) \
    && echo "sorted PAFs IDENTICAL across contig order" || echo "sorted PAFs DIFFER across contig order"
else echo "order check SKIPPED (a gold.paf is missing)"; fi
echo "historical reference: 43126"
echo "### DONE $(date -Is)"
