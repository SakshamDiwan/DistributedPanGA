#!/bin/bash
# Step 6b: unmodified FastGA 178afcb plus two diagnostic prints (private copy
# ~/pga-build/fastga-diag): the leftover bepos on entry to every same-strand
# self pair, and each thread's range of anchor contigs. Used to check the
# values the match-fastga-order build computes. Leaves gold.paf alone.
# Dataset: --export=ALL,PGA_D=<fixture dir>,PGA_STEM=<stem>,PGA_T=<threads>
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --time=01:30:00
#SBATCH --job-name=fastga-diag
#SBATCH --output=/global/u2/h/hsc53/pga-jobs/human/logs/%x.%j.log
set -uo pipefail
DIAG=$HOME/pga-build/fastga-diag/FastGA
module load cpu; module unload darshan
export PATH=$HOME/FASTGA:$PATH
RUN=$SCRATCH/pga-human/$SLURM_JOB_NAME-$SLURM_JOB_ID; mkdir -p "${RUN:?}/tmp"
export TMPDIR=$RUN/tmp
D=${PGA_D:?}; S=${PGA_STEM:?}; T=${PGA_T:-32}
echo "JOB: $SLURM_JOB_ID  NODE: $SLURM_JOB_NODELIST  DATE: $(date -Is)  RUN: $RUN"
echo "dataset: $D/$S  threads: $T  FastGA diag sha256 $(sha256sum < $DIAG | cut -c1-16)"
cd $D; [ -s $S.gix ] || { echo "index missing"; exit 1; }
/usr/bin/time -f "  FastGA %e s, %M kB maxRSS" $DIAG -1:$RUN/diag -T$T -P$TMPDIR $S 2> $RUN/diag.err || { tail -5 $RUN/diag.err; exit 1; }
tail -1 $RUN/diag.err
grep -h '^\[fga-self\]' $RUN/diag.err | LC_ALL=C sort -t= -k3 -n > $RUN/fga-self.txt
grep -h '^\[fga-range\]' $RUN/diag.err | LC_ALL=C sort -t= -k2 -n > $RUN/fga-range.txt
echo "self pairs: $(wc -l < $RUN/fga-self.txt)  ranges: $(wc -l < $RUN/fga-range.txt)"
ALNtoPAF -T$T $RUN/diag.1aln | LC_ALL=C sort -S 32G --parallel=16 -T $TMPDIR > $RUN/diag.sorted.paf
LC_ALL=C sort -S 32G --parallel=16 -T $TMPDIR gold.paf | cmp -s - $RUN/diag.sorted.paf && echo "diag output == gold" || echo "diag output DIFFERS from gold"
rm -rf "${RUN:?}/tmp" $RUN/diag.1aln $RUN/diag.sorted.paf
echo "### DONE $(date -Is)"
