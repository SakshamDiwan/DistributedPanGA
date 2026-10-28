#!/bin/bash
# Step 3 (larger datasets) with PR 1 + PR 2 (worktree ~/pga-build/pr2):
#   a. scerevisiae8 regression at N=1,2,4,8 (N=2 aborted before R-06); gold 50,046
#   b. primates4.chr6 at N=16, 8, 4 (one rep each); gold 1,547,485.
#      Only 10 contig pairs exist, so N>10 leaves ranks without alignment work.
#   c. primates14.chr6: readiness only (fixture, GDB stats), no run
# Every run: gold check (sorted, duplicates kept), [pr1] timing/memory, balance report.
# Submit: sbatch ~/pga-jobs/step3/pga_step3_l1.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --job-name=pga-step3-l1
#SBATCH --output=%x.%j.log
set -uo pipefail
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
W=$HOME/pga-build/pr2; cd "$W"; R=$W
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1 MAKEFLAGS=CC=cc
export PATH="$HOME/FASTGA:$PATH"
L=$SCRATCH/pga-step3-l1-$SLURM_JOB_ID; mkdir -p "$L"
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)  LOGS: $L"
echo "COMMIT: $(git rev-parse HEAD)  BRANCH: $(git rev-parse --abbrev-ref HEAD)  DIRTY: $(git status --porcelain | tr '\n' ' ')"
make -B CC=cc all > $L/build.log 2>&1; echo "build rc=$?  darshan: $(ldd pga-mpi-merge | grep -c darshan)"
BIN=$L/pga-mpi-merge; cp pga-mpi-merge $BIN; sha256sum $BIN | sed 's|/.*/||'
source $HOME/pga-jobs/pr1/pr1_lib.sh
F=$SCRATCH/pga-fixtures
SG=$L/sc8.gold.sorted.paf; LC_ALL=C sort $F/scerevisiae8/gold.paf > $SG
PG=$L/p4.gold.sorted.paf;  LC_ALL=C sort $F/primates4-chr6/gold.paf > $PG
echo "gold: scerevisiae8 $(wc -l < $SG)  primates4.chr6 $(wc -l < $PG)"
sha256sum $F/scerevisiae8/gold.paf $F/primates4-chr6/gold.paf | sed 's|/.*/||'

echo; echo "### a. scerevisiae8 regression ###"
for N in 2 1 4 8; do
  echo "=== sc8 N=$N $(date -Is) ==="; run_pr1 sc8 $BIN $F/scerevisiae8/scerevisiae8 $SG $N $L/sc8-N$N
done

echo; echo "### b. primates4.chr6 ###"
for N in 16 8 4; do
  echo "=== p4 N=$N $(date -Is) ==="
  ( export RUN_LIMIT_MIN=90; run_pr1 p4 $BIN $F/primates4-chr6/primates4_chr6 $PG $N $L/p4-N$N )
  grep -hoE 'stage2_rounds=[0-9]+|redist_n_rounds=[0-9]+|redist_staging=[0-9]+' $L/p4-N$N/merge.err | sort | uniq -c | head -4 | sed 's/^/  /'
done

echo; echo "### c. primates14.chr6 readiness ###"
P14=$F/primates14-chr6
ls -la $P14/primates14_chr6.1gdb $P14/.primates14_chr6.bps $P14/gold.paf 2>&1 | sed 's/^/  /'
echo "  gold alignments: $(wc -l < $P14/gold.paf)"
echo "### DONE $(date -Is)"
