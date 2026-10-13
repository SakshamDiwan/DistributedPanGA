#!/bin/bash
# PR 1 validation (branch pr1/align-instrumentation):
#   1. build; unit tier already run on the login node
#   2. small fixture: e2e at N=1,2,4,8; ASan driver at N=1,2,4 (no ASan errors, gold match)
#   3. yeast7 and scerevisiae8: main (fd4b6e9) vs PR 1 at N=1,8 (memory, overhead),
#      PR 1 at N=1,2,4,8 x3 reps with the alignment-balance report
# Submit: sbatch ~/pga-jobs/pr1/pga_pr1.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --time=01:15:00
#SBATCH --job-name=pga-pr1
#SBATCH --output=%x.%j.log
set -uo pipefail
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
R=$HOME/DistributedPanGA; cd "$R"
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1 MAKEFLAGS=CC=cc
export PATH="$HOME/FASTGA:$PATH"
L=$SCRATCH/pga-pr1-$SLURM_JOB_ID; mkdir -p "$L"
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)  LOGS: $L"
echo "COMMIT: $(git rev-parse HEAD)  BRANCH: $(git rev-parse --abbrev-ref HEAD)  DIRTY: $(git status --porcelain | tr '\n' ' ')"
make -B CC=cc all > $L/build.log 2>&1; echo "build rc=$?  darshan: $(ldd pga-mpi-merge | grep -c darshan)"
make -B CC=cc pga-mpi-merge-asan > $L/build-asan.log 2>&1; echo "asan build rc=$?"
PR1=$L/pga-mpi-merge.pr1; cp pga-mpi-merge $PR1
ASAN=$L/pga-mpi-merge.asan; cp tests/bin/pga-mpi-merge-asan $ASAN
BASE=$HOME/pga-build/main-fd4b6e9/pga-mpi-merge
sha256sum $BASE $PR1 $ASAN | sed 's|/.*/||'
source $HOME/pga-jobs/pr1/pr1_lib.sh

Y=$R/tests/fixtures/yeast;               YG=$L/yeast7.gold.sorted.paf; LC_ALL=C sort $Y/gold.paf > $YG
S=$SCRATCH/pga-fixtures/scerevisiae8;     SG=$L/sc8.gold.sorted.paf;   LC_ALL=C sort $S/gold.paf > $SG
SM=$R/tests/fixtures/generated;           MG=$L/small.gold.sorted.paf; LC_ALL=C sort $SM/gold.paf > $MG
echo "gold: yeast7 $(wc -l < $YG)  scerevisiae8 $(wc -l < $SG)  small6 $(wc -l < $MG)"
sha256sum $Y/gold.paf $S/gold.paf | sed 's|/.*/||'

echo; echo "### 2. small fixture ###"
RANK_COUNTS="1 2 4 8" bash scripts/verify_e2e_small.sh > $L/e2e-small.log 2>&1; echo "e2e-small rc=$?"; grep -E 'PASS|FAIL' $L/e2e-small.log | tail -1
for N in 1 2 4; do
  ( export ASAN_OPTIONS=detect_leaks=0; run_pr1 asan $ASAN $SM/small6 $MG $N $L/asan-small-N$N ) | grep -E '^RESULT'
  echo "  ASan errors: $(grep -c 'ERROR: AddressSanitizer' $L/asan-small-N$N/merge.err)   UBSan outside lib/: $(grep 'runtime error' $L/asan-small-N$N/merge.err | grep -vc '^lib/')"
done

for DS in yeast7 sc8; do
  if [ $DS = yeast7 ]; then STM=$Y/yeast7; G=$YG; else STM=$S/scerevisiae8; G=$SG; fi
  echo; echo "### 3. $DS ###"
  for rep in 1 2 3; do
    for N in 1 8; do
      echo "=== $DS rep=$rep N=$N base $(date -Is) ==="; run_pr1 base $BASE $STM $G $N $L/$DS-base-N$N-r$rep
    done
    for N in 1 2 4 8; do
      echo "=== $DS rep=$rep N=$N pr1 $(date -Is) ==="; run_pr1 pr1 $PR1 $STM $G $N $L/$DS-pr1-N$N-r$rep
    done
  done
done
echo "### DONE $(date -Is)"
