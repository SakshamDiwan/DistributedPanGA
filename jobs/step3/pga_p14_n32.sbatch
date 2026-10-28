#!/bin/bash
# primates14.chr6 at N=32 on 2 nodes (16 ranks/node; binary addf993); gold 4,577,556.
# First multi-node run of the pipeline. Scratch files on Lustre.
# Submit: sbatch ~/pga-jobs/step3/pga_p14_n32.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --job-name=pga-p14-n32
#SBATCH --output=%x.%j.log
set -uo pipefail
L=$SCRATCH/pga-p14-n32-$SLURM_JOB_ID; mkdir -p "$L"
source $HOME/pga-jobs/step3/diag_lib.sh; header
P14G=$L/p14.gold.sorted.paf; LC_ALL=C sort -S 8G $F/primates14-chr6/gold.paf > $P14G
echo "gold: primates14 $(wc -l < $P14G)"
mkdir -p $L/scratch
echo "=== p14 N=32 $(date -Is) ==="
( export RUN_LIMIT_MIN=160 PGA_SCRATCH_DIR=$L/scratch; run_pr1 p14 $BIN $F/primates14-chr6/primates14_chr6 $P14G 32 $L/p14-N32 )
grep -hoE 'stage2_rounds=[0-9]+|redist_n_rounds=[0-9]+' $L/p14-N32/merge.err | sort | uniq -c | head -3 | sed 's/^/  /'
tail -1 $L/p14-N32/pga-mpi-merge.timing.csv | sed 's/^/  csv: /'
sed -n '/== 4/,$p' $L/p14-N32/balance.txt | head -40 | sed 's/^/  /'
echo "### DONE $(date -Is)"
