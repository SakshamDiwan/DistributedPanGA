#!/bin/bash
# Week 1, Job D: extraction baseline on yeast7. 1 node, N in {1,2,4,8}, 3 reps, interleaved.
# Fixed per-rank resources: -c 16 (8 physical cores), --cpu-bind=cores, OMP_NUM_THREADS=1,
# FREQ=10, merge NTHREADS=1, sort nthreads=4 (compiled in). Every run is also checked vs gold.
# Submit ONLY after Job C passes: sbatch ~/pga-jobs/week1/pga_w1_D.sbatch   (YEAST=yeast-alpha via --export if needed)
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --job-name=pga-w1-D
#SBATCH --output=%x.%j.log
set -uo pipefail
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
R=$HOME/DistributedPanGA; cd "$R"
YEAST=${YEAST:-yeast}
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1
export PATH="$HOME/FASTGA:$PATH"
STEM=$R/tests/fixtures/$YEAST/yeast7
L=$SCRATCH/pga-w1-D-$SLURM_JOB_ID; mkdir -p "$L"
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)  YEAST=$YEAST  LOGS: $L"
echo "COMMIT: $(git rev-parse HEAD)  DIRTY: $(git status --porcelain | tr '\n' ' ')"
module -t list 2>&1 | tr '\n' ' '; echo
echo "BUILD: $(make -Bn CC=cc pga-mpi-merge | head -1 | cut -c1-80)..."
sha256sum pga-mpi-merge $STEM.fa.gz $STEM.1gdb ${STEM%/yeast7}/gold.paf
lscpu | grep -E 'Model name|^CPU\(s\)|Thread|Socket|NUMA node\(s\)'
GOLD_SORTED=$L/gold.sorted.paf; LC_ALL=C sort ${STEM%/yeast7}/gold.paf > $GOLD_SORTED
source $HOME/pga-jobs/week1/pga_w1_e2e_lib.sh
for rep in 1 2 3; do
  for N in 1 2 4 8; do
    echo; echo "=== rep=$rep N=$N $(date -Is) ==="
    run_one $N $L/N$N-r$rep
  done
done
echo; echo "### summary"; grep -h '^RESULT' $HOME/pga-jobs/week1/logs/pga-w1-D.$SLURM_JOB_ID.log 2>/dev/null | wc -l
echo; sacct -j $SLURM_JOB_ID --format=JobID%-16,JobName%-14,Elapsed,MaxRSS,AveRSS,NTasks,State -P
echo "### DONE $(date -Is)"
