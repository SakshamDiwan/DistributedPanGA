#!/bin/bash
# O1a rerun (A/B only). Gates test-mpi, phase1/2, stage1/2 already passed on
# f528988 in job 59686567, which then hung on its third timed run (base binary,
# GDB on /tmp, N=1; >32 min vs 1:53 on Lustre). Changes:
#   - every timed srun has its own limit (RUN_LIMIT_MIN, default 10 min);
#   - Lustre runs (the main A/B) come first, the /tmp control runs after;
#   - LCP gate stays last.
#   base = pre-O1a pga-mpi-merge (abaea3f, worktree ~/pga-build/o1-0-abaea3f)
#   o1a  = this branch, rebuilt here
#   H2 control: GDB on Lustre ($SCRATCH) vs node-local /tmp (tmpfs), N in {1,8}
#   plus o1a at N in {2,4} on Lustre to complete the scaling curve. 3 reps, interleaved.
# Submit: sbatch ~/pga-jobs/o1/pga_o1_a2.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --time=01:30:00
#SBATCH --job-name=pga-o1-a2
#SBATCH --output=%x.%j.log
set -uo pipefail
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
R=$HOME/DistributedPanGA; cd "$R"
BASE_BIN=$HOME/pga-build/o1-0-abaea3f/pga-mpi-merge
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1 MAKEFLAGS=CC=cc
export PATH="$HOME/FASTGA:$PATH"
export PGA_FIXTURES=$R/tests/fixtures/yeast PGA_GDB_STEM=$R/tests/fixtures/yeast/yeast7
LUSTRE_STEM=$PGA_GDB_STEM
L=$SCRATCH/pga-o1-a2-$SLURM_JOB_ID; mkdir -p "$L"
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)  LOGS: $L"
echo "COMMIT: $(git rev-parse HEAD)  BRANCH: $(git rev-parse --abbrev-ref HEAD)  DIRTY: $(git status --porcelain | tr '\n' ' ')"
echo "BASE:   $(git -C $HOME/pga-build/o1-0-abaea3f rev-parse HEAD)"
module -t list 2>&1 | tr '\n' ' '; echo
make -B CC=cc all > $L/build.log 2>&1; echo "build rc=$?  darshan-linked: $(ldd pga-mpi-merge | grep -c darshan)"
echo "BUILD: $(make -Bn CC=cc pga-mpi-merge | head -1 | cut -c1-100)..."
O1A_BIN=$L/pga-mpi-merge.o1a; cp pga-mpi-merge $O1A_BIN    # frozen copy: the LCP verifier runs make clean
sha256sum $BASE_BIN $O1A_BIN $LUSTRE_STEM.fa.gz $LUSTRE_STEM.1gdb $PGA_FIXTURES/gold.paf
GOLD_SORTED=$L/gold.sorted.paf; LC_ALL=C sort $PGA_FIXTURES/gold.paf > $GOLD_SORTED
source $HOME/pga-jobs/o1/o1_0_summary.sh
source $HOME/pga-jobs/o1/o1_lib.sh

TMPD=/tmp/pga-gdb-$SLURM_JOB_ID; [ "${RUN_TMP_CONTROL:-0}" = 1 ] && mkdir -p $TMPD
[ "${RUN_TMP_CONTROL:-0}" = 1 ] && cp $LUSTRE_STEM.1gdb $(dirname $LUSTRE_STEM)/.yeast7.bps $TMPD/
TMP_STEM=$TMPD/yeast7
echo "tmp copy: $(df -T $TMPD | tail -1 | awk '{print $2}')  $(cmp -s $TMPD/.yeast7.bps $(dirname $LUSTRE_STEM)/.yeast7.bps && echo bps-identical)"

echo; echo "### 1. gates: test-mpi, phase1/2, stage1/2 passed on this commit in job 59686567 (not rerun) ###"

echo; echo "### 2a. A/B on Lustre: base vs o1a at N=1,8; o1a at N=2,4; 3 reps interleaved ###"
for rep in 1 2 3; do
  for N in 1 8; do for b in base o1a; do
    [ $b = base ] && BIN=$BASE_BIN || BIN=$O1A_BIN
    echo "=== rep=$rep N=$N fs=lustre bin=$b $(date -Is) ==="
    run_bin "$b/lustre" $BIN $LUSTRE_STEM $N $L/$b-lustre-N$N-r$rep
  done; done
  for N in 2 4; do
    echo "=== rep=$rep N=$N fs=lustre bin=o1a $(date -Is) ==="
    run_bin "o1a/lustre" $O1A_BIN $LUSTRE_STEM $N $L/o1a-lustre-N$N-r$rep
  done
done

echo; echo "### 2b. H2 control: GDB on node-local /tmp, N=1,8, 3 reps (o1a first) ###"
# Off by default: job 59686567 hung (node left in CG) on its first /tmp run, cause
# unknown. Enable only deliberately: sbatch --export=ALL,RUN_TMP_CONTROL=1 ...
if [ "${RUN_TMP_CONTROL:-0}" != 1 ]; then echo "skipped (RUN_TMP_CONTROL!=1)"; fi
for rep in $([ "${RUN_TMP_CONTROL:-0}" = 1 ] && echo 1 2 3); do
  for N in 8 1; do for b in o1a base; do
    [ $b = base ] && BIN=$BASE_BIN || BIN=$O1A_BIN
    echo "=== rep=$rep N=$N fs=tmp bin=$b $(date -Is) ==="
    run_bin "$b/tmp" $BIN $TMP_STEM $N $L/$b-tmp-N$N-r$rep
  done; done
done

echo; echo "### 3. verify-phase3-stage2-lcp (runs make clean; last) ###"
t0=$(date +%s); RANK_COUNTS="1 2 4 8" make verify-phase3-stage2-lcp > $L/lcp.log 2>&1; echo "lcp rc=$? ($(( $(date +%s)-t0 )) s)"
grep -E 'ALL PASS|FAIL' $L/lcp.log | tail -3
make -B CC=cc all > /dev/null 2>&1; echo "rebuild rc=$?  same-as-measured: $(cmp -s pga-mpi-merge $O1A_BIN && echo yes || echo NO)"
rm -rf $TMPD
echo "### DONE $(date -Is)"
