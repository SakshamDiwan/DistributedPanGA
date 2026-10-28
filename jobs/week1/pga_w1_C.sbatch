#!/bin/bash
# Week 1, Job C: yeast7 phase verifiers + end-to-end vs independent FastGA gold.
# Candidate dir: YEAST=yeast (default, file order) or YEAST=yeast-alpha.
# Submit: sbatch ~/pga-jobs/week1/pga_w1_C.sbatch      (or: sbatch --export=ALL,YEAST=yeast-alpha ~/pga-jobs/week1/pga_w1_C.sbatch)
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --job-name=pga-w1-C
#SBATCH --output=%x.%j.log
set -uo pipefail
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
R=$HOME/DistributedPanGA; cd "$R"
YEAST=${YEAST:-yeast}
export MPICH_GPU_SUPPORT_ENABLED=0 MAKEFLAGS=CC=cc OMP_NUM_THREADS=1
export PATH="$HOME/FASTGA:$PATH"
export PGA_FIXTURES="$R/tests/fixtures/$YEAST"
export PGA_GDB_STEM="$PGA_FIXTURES/yeast7"
STEM=$PGA_GDB_STEM
L=$SCRATCH/pga-w1-C-$SLURM_JOB_ID; mkdir -p "$L"
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)  YEAST=$YEAST  LOGS: $L"
echo "COMMIT: $(git rev-parse HEAD)  DIRTY: $(git status --porcelain | tr '\n' ' ')"
module -t list 2>&1 | tr '\n' ' '; echo
sha256sum pga-extract pga-sort pga-mpi-stage1 pga-mpi pga-mpi-merge
( cd $PGA_FIXTURES && sha256sum yeast7.fa.gz yeast7.1gdb gold.paf gold.tuples 2>&1 )

echo; echo "### 1. phase verifiers on yeast7 ###"
for t in verify-phase1 verify-phase2; do
  t0=$(date +%s); make $t > $L/$t.log 2>&1; echo "$t rc=$? ($(( $(date +%s)-t0 )) s)"; grep -E 'PASS|FAIL' $L/$t.log | tail -3
done
for t in verify-phase3-stage1 verify-phase3-stage2; do
  t0=$(date +%s); RANK_COUNTS="1 2 4 8" make $t > $L/$t.log 2>&1; echo "$t rc=$? ($(( $(date +%s)-t0 )) s)"; grep -E 'PASS|FAIL' $L/$t.log | tail -9
done

echo; echo "### 2. end-to-end: full sorted aggregate PAF vs gold, 1 2 4 8 ranks ###"
GOLD_SORTED=$L/gold.sorted.paf; LC_ALL=C sort $PGA_FIXTURES/gold.paf > $GOLD_SORTED
echo "gold: $(wc -l < $GOLD_SORTED) alignments, $(uniq -d $GOLD_SORTED | wc -l) duplicated lines  (historical: 43126)"
source $HOME/pga-jobs/week1/pga_w1_e2e_lib.sh
for N in 1 2 4 8; do run_one $N $SCRATCH/pga-w1-C-$SLURM_JOB_ID/e2e-N$N; done
echo "consistency:"; for N in 2 4 8; do cmp -s $L/e2e-N1/got.sorted.paf $L/e2e-N$N/got.sorted.paf && echo "  N=$N identical to N=1" || echo "  N=$N DIFFERS from N=1"; done

echo; echo "### 3. verify-phase3-stage2-lcp (runs make clean; last) ###"
t0=$(date +%s); RANK_COUNTS="1 2 4 8" make verify-phase3-stage2-lcp > $L/lcp.log 2>&1; echo "lcp rc=$? ($(( $(date +%s)-t0 )) s)"
grep -E 'PASS|FAIL' $L/lcp.log | tail -9
make -B all > /dev/null 2>&1; echo "rebuild rc=$?"; sha256sum pga-extract pga-sort pga-mpi-stage1 pga-mpi pga-mpi-merge
echo; sacct -j $SLURM_JOB_ID --format=JobID%-16,JobName%-14,Elapsed,MaxRSS,NTasks,State -P | head -40
echo "### DONE $(date -Is)"
