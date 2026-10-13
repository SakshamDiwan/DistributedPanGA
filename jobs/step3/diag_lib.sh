# Shared by the diagnostic / primates14 jobs. Uses the prebuilt diag binary.
BIN=$HOME/pga-jobs/step3/bin/pga-mpi-merge.addf993
R=$HOME/pga-build/pr2        # for scripts/analyze_align_balance.py
F=$SCRATCH/pga-fixtures
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1 PATH="$HOME/FASTGA:$PATH"
source $HOME/pga-jobs/pr1/pr1_lib.sh
header() {
  echo "JOB: $SLURM_JOB_ID  NODES: $SLURM_JOB_NODELIST  DATE: $(date -Is)  LOGS: $L"
  echo "BIN: $(sha256sum $BIN | cut -c1-16)  source commit addf993 (diag/align-cost)"
  echo "/tmp: $(df -h /tmp | tail -1 | awk '{print $1, $2, "size", $4, "free"}')"
}
