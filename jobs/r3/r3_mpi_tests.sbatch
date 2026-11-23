#!/bin/bash
# R3: MPI test tier (tests/run_mpi.sh) on a compute node, incl. M-07 (run_stage2_take, sort settings)
# and the MPI_Alltoallv_c cases of test_chunked_mpi. Work dir on $SCRATCH (visible to every rank).
# Fails if any M-06 digest is missing, so the cross-world-size check cannot pass vacuously.
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:15:00
#SBATCH --job-name=r3mpitest
#SBATCH --output=/global/u2/h/hsc53/pga-jobs/r3/logs/%x.%j.log
set -uo pipefail
module load cpu; module unload darshan
export MPICH_GPU_SUPPORT_ENABLED=0
R=$HOME/pga-build/r3
echo "repo: $(git -C $R rev-parse HEAD)  node: $SLURM_JOB_NODELIST"
export TMPDIR=$SCRATCH/pga-r3/mpitest-$SLURM_JOB_ID; mkdir -p "$TMPDIR"
export RANK_COUNTS="1 2 3 4 8"
bash $R/tests/run_mpi.sh 2>&1 | tee $TMPDIR/out.txt
rc=${PIPESTATUS[0]}
m=$(grep -c 'matches ws=' $TMPDIR/out.txt)
echo "M-06 comparisons made: $m (want 4)"
[ "$rc" = 0 ] && [ "$m" = 4 ] && echo "MPI TIER: PASS" || echo "MPI TIER: FAIL (rc=$rc, M-06 comparisons $m)"
rm -rf "${TMPDIR:?}"
