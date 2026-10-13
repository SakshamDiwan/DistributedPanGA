#!/bin/bash
# PR 2 validation (worktree ~/pga-build/pr2, branch pr2/large-count-exchange, stacked on PR 1):
#   1. test-mpi at ws 1..8, default budget and a 256-byte budget (multi-round stage 2)
#   2. small6 e2e at N=1,2,4,8: default and 4096-byte budget
#   3. yeast7 verify-phase3-stage2 (Tier 1 + 2) with a 1 MiB budget at N=1,2,4,8
#   4. yeast7 / scerevisiae8 e2e with an 8 MiB budget (stage-2 AND seed rounds) vs PR 1
#      binary at default: gold match AND raw per-rank PAF byte-identical
# Submit: sbatch ~/pga-jobs/pr2/pga_pr2.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --job-name=pga-pr2
#SBATCH --output=%x.%j.log
set -uo pipefail
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
W=$HOME/pga-build/pr2; cd "$W"
MAIN=$HOME/DistributedPanGA
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1 MAKEFLAGS=CC=cc
export PATH="$HOME/FASTGA:$PATH"
L=$SCRATCH/pga-pr2-$SLURM_JOB_ID; mkdir -p "$L"
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)  LOGS: $L"
echo "COMMIT: $(git rev-parse HEAD)  BRANCH: $(git rev-parse --abbrev-ref HEAD)  DIRTY: $(git status --porcelain | tr '\n' ' ')"
make -B CC=cc all > $L/build.log 2>&1; echo "build rc=$?  darshan: $(ldd pga-mpi-merge | grep -c darshan)"
PR2=$L/pga-mpi-merge.pr2; cp pga-mpi-merge $PR2
# PR 1 reference binary, built from its commit in a scratch worktree-free way
git -C $MAIN show pr1/align-instrumentation:src/pga-mpi-merge.c > /dev/null || echo "WARN: PR1 branch missing"
PR1=$L/pga-mpi-merge.pr1
( T=$(mktemp -d -p $L) && git -C $MAIN archive pr1/align-instrumentation | tar -x -C $T && cd $T && make CC=cc pga-mpi-merge > $L/build-pr1.log 2>&1 && cp pga-mpi-merge $PR1 ) && echo "pr1 reference built" || echo "pr1 reference build FAILED"
sha256sum $PR1 $PR2 | sed 's|/.*/||'
source $HOME/pga-jobs/pr1/pr1_lib.sh
R=$W

echo; echo "### 1. test-mpi ###"
RANK_COUNTS="1 2 3 4 8" make test-mpi > $L/test-mpi.log 2>&1; echo "test-mpi default rc=$?  FAIL lines: $(grep -c FAIL $L/test-mpi.log)"
PGA_EXCHANGE_ROUND_BYTES=256 RANK_COUNTS="1 2 3 4 8" make test-mpi > $L/test-mpi-256.log 2>&1; echo "test-mpi budget=256 rc=$?  FAIL lines: $(grep -c FAIL $L/test-mpi-256.log)"

echo; echo "### 2. small6 e2e ###"
export PGA_GDB_STEM=$MAIN/tests/fixtures/generated/small6 PGA_GOLD_PAF=$MAIN/tests/fixtures/generated/gold.paf
RANK_COUNTS="1 2 4 8" bash scripts/verify_e2e_small.sh > $L/e2e-small.log 2>&1; echo "e2e small default rc=$?"; grep -E 'ALL PASS|FAIL' $L/e2e-small.log | tail -1
PGA_EXCHANGE_ROUND_BYTES=4096 RANK_COUNTS="1 2 4 8" bash scripts/verify_e2e_small.sh > $L/e2e-small-4k.log 2>&1; echo "e2e small budget=4096 rc=$?"; grep -E 'ALL PASS|FAIL' $L/e2e-small-4k.log | tail -1
unset PGA_GDB_STEM PGA_GOLD_PAF

echo; echo "### 3. yeast7 verify-phase3-stage2, 1 MiB budget ###"
t0=$(date +%s)
PGA_FIXTURES=$MAIN/tests/fixtures/yeast PGA_GDB_STEM=$MAIN/tests/fixtures/yeast/yeast7 PGA_EXCHANGE_ROUND_BYTES=1048576 \
  RANK_COUNTS="1 2 4 8" make verify-phase3-stage2 > $L/stage2-1m.log 2>&1
echo "verify-phase3-stage2 budget=1MiB rc=$? ($(( $(date +%s)-t0 )) s)"; grep -E 'ALL PASS|FAIL' $L/stage2-1m.log | tail -3

echo; echo "### 4. e2e: PR2 at 8 MiB budget vs PR1 default ###"
Y=$MAIN/tests/fixtures/yeast; S=$SCRATCH/pga-fixtures/scerevisiae8
YG=$L/yeast7.gold.sorted.paf; LC_ALL=C sort $Y/gold.paf > $YG
SG=$L/sc8.gold.sorted.paf;   LC_ALL=C sort $S/gold.paf > $SG
for spec in "yeast7 $Y/yeast7 $YG 1" "yeast7 $Y/yeast7 $YG 2" "yeast7 $Y/yeast7 $YG 8" "sc8 $S/scerevisiae8 $SG 4"; do
  set -- $spec
  echo "=== $1 N=$4 ==="
  run_pr1 pr1 $PR1 $2 $3 $4 $L/$1-pr1-N$4
  ( export PGA_EXCHANGE_ROUND_BYTES=8388608; run_pr1 pr2-8M $PR2 $2 $3 $4 $L/$1-pr2-N$4 )
  grep -hoE 'stage2_rounds=[0-9]+|redist_n_rounds=[0-9]+|redist_c_rounds=[0-9]+|redist_staging=[0-9]+' $L/$1-pr2-N$4/merge.err | sort | uniq -c | sed 's/^/  /' | head -8
  same=1; for f in $L/$1-pr1-N$4/out.rank*.paf; do cmp -s "$f" "$L/$1-pr2-N$4/$(basename $f)" || same=0; done
  echo "  raw per-rank PAF PR2(8 MiB) vs PR1(default): $([ $same = 1 ] && echo IDENTICAL || echo DIFFERENT)"
done
echo "### DONE $(date -Is)"
