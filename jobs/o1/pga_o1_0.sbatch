#!/bin/bash
# O1-0: instrumented baseline + pre-sort shuffle experiment on yeast7.
#   1. correctness for this build: test-mpi, verify-phase3-stage2 (pga-mpi links the changed dsort.c)
#   2. instrumented baseline: N in {1,2,4,8} x 3 reps, Job D configuration, gold check each run
#   3. shuffle: PGA_PRESORT_SHUFFLE in {1,2,3} at N=1 and N=4, gold check + raw per-rank PAF cmp vs unshuffled
# Submit: sbatch ~/pga-jobs/o1/pga_o1_0.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --time=01:15:00
#SBATCH --job-name=pga-o1-0
#SBATCH --output=%x.%j.log
set -uo pipefail
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
R=$HOME/DistributedPanGA; cd "$R"
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1 MAKEFLAGS=CC=cc
export PATH="$HOME/FASTGA:$PATH"
export PGA_FIXTURES=$R/tests/fixtures/yeast PGA_GDB_STEM=$R/tests/fixtures/yeast/yeast7
STEM=$PGA_GDB_STEM
L=$SCRATCH/pga-o1-0-$SLURM_JOB_ID; mkdir -p "$L"
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)  LOGS: $L"
echo "COMMIT: $(git rev-parse HEAD)  BRANCH: $(git rev-parse --abbrev-ref HEAD)  DIRTY: $(git status --porcelain | tr '\n' ' ')"
module -t list 2>&1 | tr '\n' ' '; echo
echo "BUILD: $(make -Bn CC=cc pga-mpi-merge | head -1 | cut -c1-100)..."
sha256sum pga-mpi pga-mpi-merge $STEM.fa.gz $STEM.1gdb $PGA_FIXTURES/gold.paf
GOLD_SORTED=$L/gold.sorted.paf; LC_ALL=C sort $PGA_FIXTURES/gold.paf > $GOLD_SORTED
source $HOME/pga-jobs/week1/pga_w1_e2e_lib.sh
source $HOME/pga-jobs/o1/o1_0_summary.sh

echo; echo "### 1. correctness for this build ###"
RANK_COUNTS="1 2 3 4 8" make test-mpi > $L/test-mpi.log 2>&1; echo "test-mpi rc=$?"; grep -cE 'FAIL' $L/test-mpi.log | sed 's/^/  FAIL lines: /'
t0=$(date +%s); RANK_COUNTS="1 2 4 8" make verify-phase3-stage2 > $L/stage2.log 2>&1
echo "verify-phase3-stage2 rc=$? ($(( $(date +%s)-t0 )) s)"; grep -E 'PASS|FAIL' $L/stage2.log | tail -9

echo; echo "### 2. instrumented baseline (Job D config) ###"
for rep in 1 2 3; do for N in 1 2 4 8; do
  echo "=== rep=$rep N=$N $(date -Is) ==="
  run_one $N $L/base-N$N-r$rep
  o1_summary $L/base-N$N-r$rep
done; done

echo; echo "### 3. pre-sort shuffle experiment ###"
for N in 1 4; do for seed in 1 2 3; do
  echo "=== shuffle seed=$seed N=$N $(date -Is) ==="
  export PGA_PRESORT_SHUFFLE=$seed
  run_one $N $L/shuf-N$N-s$seed
  unset PGA_PRESORT_SHUFFLE
  grep -c 'PGA_PRESORT_SHUFFLE' $L/shuf-N$N-s$seed/merge.err | sed 's/^/  shuffle banner lines: /'
  o1_summary $L/shuf-N$N-s$seed
  same=1; for f in $L/base-N$N-r1/out.rank*.paf; do cmp -s "$f" "$L/shuf-N$N-s$seed/$(basename $f)" || same=0; done
  echo "  raw per-rank PAF vs unshuffled rep1: $([ $same = 1 ] && echo IDENTICAL || echo DIFFERENT)"
done; done

echo; sacct -j $SLURM_JOB_ID --format=JobID%-16,JobName%-20,Elapsed,State -P | tail -n +2 | grep -vc COMPLETED | sed 's/^/non-COMPLETED steps (incl. running batch): /'
echo "### DONE $(date -Is)"
