#!/bin/bash
# Week 1, Job A: MPI tier + small-fixture phase verifiers + e2e vs FastGA gold.
# Verifies merged main (b4f5d3b) built with cc. Submit: sbatch ~/pga-jobs/week1/pga_w1_A.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=32
#SBATCH --time=00:25:00
#SBATCH --job-name=pga-w1-A
#SBATCH --output=%x.%j.log
set -uo pipefail
# cc links darshan by default; darshan leaves errno=95 set, which trips ONElib's
# post-fclose errno check (lib/ONElib.c:412) in every MPI binary. CPU-only env.
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
module -t list 2>&1 | tr "\n" " "; echo
R=$HOME/DistributedPanGA; cd "$R"
export MPICH_GPU_SUPPORT_ENABLED=0
export MAKEFLAGS=CC=cc                 # keep every (re)build on cc, incl. the LCP verifier's make clean
export FASTGA_BIN="$HOME/FASTGA"; export PATH="$FASTGA_BIN:$PATH"
export PGA_FIXTURES="$R/tests/fixtures/generated"
export PGA_GDB_STEM="$PGA_FIXTURES/small6"
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)"
echo "COMMIT: $(git rev-parse HEAD)  BRANCH: $(git rev-parse --abbrev-ref HEAD)"
echo "DIRTY:  $(git status --porcelain | tr "\n" " ")"
L=$SCRATCH/pga-w1-A-$SLURM_JOB_ID; mkdir -p "$L"; echo "LOGS: $L"
echo "### binary hashes (start)"; sha256sum pga-extract pga-sort pga-mpi-stage1 pga-mpi pga-merge pga-mpi-merge

echo; echo "### 1. MPI tier (1 2 3 4 8) ###"
RANK_COUNTS="1 2 3 4 8" make test-mpi 2>&1 | tail -60
echo "### test-mpi rc=${PIPESTATUS[0]}"

echo; echo "### 2. phase verifiers, small fixture ###"
for t in verify-phase1 verify-phase2; do
  make $t > $L/$t.log 2>&1; echo "$t rc=$?"; grep -E 'PASS|FAIL' $L/$t.log | tail -3
done
for t in verify-phase3-stage1 verify-phase3-stage2; do
  RANK_COUNTS="1 2 4 8" make $t > $L/$t.log 2>&1; echo "$t rc=$?"; grep -E 'PASS|FAIL' $L/$t.log | tail -6
done

echo; echo "### 3. e2e vs FastGA gold (1 2 4 8) ###"
RANK_COUNTS="1 2 4 8" bash scripts/verify_e2e_small.sh > $L/e2e.log 2>&1; echo "### e2e rc=$?"; tail -30 $L/e2e.log

echo; echo "### 4. verify-phase3-stage2-lcp (runs make clean; last) ###"
RANK_COUNTS="1 2 4 8" make verify-phase3-stage2-lcp > $L/lcp.log 2>&1; echo "lcp rc=$?"
grep -E 'PASS|FAIL' $L/lcp.log | tail -8

echo; echo "### restore full cc build"; make -B all > /dev/null 2>&1; echo "rebuild rc=$?"
echo "### binary hashes (end)"; sha256sum pga-extract pga-sort pga-mpi-stage1 pga-mpi pga-merge pga-mpi-merge
echo "### DONE $(date -Is)"
