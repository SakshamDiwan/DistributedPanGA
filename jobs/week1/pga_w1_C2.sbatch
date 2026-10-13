#!/bin/bash
# Week 1, Job C2: rerun verify-phase3-stage2 on yeast7 after the Tier-2 awk fix
# (hex k-mers like "...83e9" were compared numerically; now forced to string compare).
# Submit: sbatch ~/pga-jobs/week1/pga_w1_C2.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=32
#SBATCH --time=00:45:00
#SBATCH --job-name=pga-w1-C2
#SBATCH --output=%x.%j.log
set -uo pipefail
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
R=$HOME/DistributedPanGA; cd "$R"
YEAST=${YEAST:-yeast}
export MPICH_GPU_SUPPORT_ENABLED=0 MAKEFLAGS=CC=cc OMP_NUM_THREADS=1
export PATH="$HOME/FASTGA:$PATH"
export PGA_FIXTURES="$R/tests/fixtures/$YEAST"
export PGA_GDB_STEM="$PGA_FIXTURES/yeast7"
L=$SCRATCH/pga-w1-C2-$SLURM_JOB_ID; mkdir -p "$L"
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)  YEAST=$YEAST  LOGS: $L"
echo "COMMIT: $(git rev-parse HEAD)  DIRTY: $(git status --porcelain | tr '\n' ' ')"
git diff scripts/verify_phase3_stage2.sh
sha256sum pga-mpi
( cd $PGA_FIXTURES && sha256sum yeast7.fa.gz gold.tuples )

t0=$(date +%s); RANK_COUNTS="1 2 4 8" make verify-phase3-stage2 > $L/verify-phase3-stage2.log 2>&1
echo "verify-phase3-stage2 rc=$? ($(( $(date +%s)-t0 )) s)"
grep -E '===|PASS|FAIL' $L/verify-phase3-stage2.log
sha256sum pga-mpi
echo "### DONE $(date -Is)"
