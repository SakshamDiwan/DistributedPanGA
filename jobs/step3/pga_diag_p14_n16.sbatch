#!/bin/bash
# Diagnostic + primates14.chr6 on 1 node (binary addf993, diag/align-cost):
#   1. primates4.chr6 at N=8 and N=16 with the new per-task features (gold 1,547,485)
#   2. primates14.chr6 at N=16 (gold 4,577,556), scratch files on Lustre
# Submit: sbatch ~/pga-jobs/step3/pga_diag_p14_n16.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --job-name=pga-diag-p14-n16
#SBATCH --output=%x.%j.log
set -uo pipefail
L=$SCRATCH/pga-diag-p14-n16-$SLURM_JOB_ID; mkdir -p "$L"
source $HOME/pga-jobs/step3/diag_lib.sh; header
PG=$L/p4.gold.sorted.paf;  LC_ALL=C sort $F/primates4-chr6/gold.paf > $PG
P14G=$L/p14.gold.sorted.paf; LC_ALL=C sort -S 8G $F/primates14-chr6/gold.paf > $P14G
echo "gold: primates4 $(wc -l < $PG)  primates14 $(wc -l < $P14G)"

echo; echo "### 1. primates4.chr6 diagnostic ###"
for N in 16 8; do
  echo "=== p4 N=$N $(date -Is) ==="
  ( export RUN_LIMIT_MIN=60; run_pr1 p4 $BIN $F/primates4-chr6/primates4_chr6 $PG $N $L/p4-N$N )
  sed -n '/== 4/,$p' $L/p4-N$N/balance.txt | head -40 | sed 's/^/  /'
done

echo; echo "### 2. primates14.chr6 N=16 ###"
mkdir -p $L/scratch
echo "=== p14 N=16 $(date -Is) ==="
( export RUN_LIMIT_MIN=200 PGA_SCRATCH_DIR=$L/scratch; run_pr1 p14 $BIN $F/primates14-chr6/primates14_chr6 $P14G 16 $L/p14-N16 )
grep -hoE 'stage2_rounds=[0-9]+|redist_n_rounds=[0-9]+' $L/p14-N16/merge.err | sort | uniq -c | head -3 | sed 's/^/  /'
tail -1 $L/p14-N16/pga-mpi-merge.timing.csv | sed 's/^/  csv: /'
sed -n '/== 4/,$p' $L/p14-N16/balance.txt | head -40 | sed 's/^/  /'
echo "### DONE $(date -Is)"
