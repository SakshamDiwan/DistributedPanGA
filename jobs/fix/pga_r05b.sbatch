#!/bin/bash
# R-05 recheck of the failure cases. In job 59694749 `srun --export=ALL,PATH=...`
# did not change the PATH the ranks used (env var PGA_ALNTOPAF_TIMEOUT did arrive),
# so the real ALNtoPAF ran. This job first shows what PATH a rank actually sees,
# then reruns the cases with PATH changed in the launching shell instead.
# Submit: sbatch ~/pga-jobs/fix/pga_r05b.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=32
#SBATCH --time=00:10:00
#SBATCH --job-name=pga-r05b
#SBATCH --output=%x.%j.log
set -uo pipefail
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
R=$HOME/DistributedPanGA; cd "$R"
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1
L=$SCRATCH/pga-r05b-$SLURM_JOB_ID; mkdir -p "$L"
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)  LOGS: $L"
echo "COMMIT: $(git rev-parse HEAD)  BRANCH: $(git rev-parse --abbrev-ref HEAD)"
BIN=$L/pga-mpi-merge; cp pga-mpi-merge $BIN; sha256sum $BIN
SMALL=$R/tests/fixtures/generated/small6
SRUN=$(command -v srun); echo "srun: $SRUN"
FAKE=$L/fakebin; mkdir -p $FAKE
printf '#!/bin/sh\necho "fake ALNtoPAF: hanging on purpose" >&2\nsleep 600\n' > $FAKE/ALNtoPAF; chmod +x $FAKE/ALNtoPAF

echo; echo "### 0. what PATH / ALNtoPAF does a rank see? ###"
echo "batch shell: $(command -v ALNtoPAF || echo none)"
echo "A (--export=ALL,PATH=...):"
$SRUN -n 1 --export=ALL,PATH=$FAKE:/usr/bin:/bin sh -c 'echo "  PATH=$PATH"; echo "  ALNtoPAF -> $(command -v ALNtoPAF || echo none)"'
echo "B (PATH set in launching shell):"
( export PATH=$FAKE:/usr/bin:/bin; $SRUN -n 1 sh -c 'echo "  PATH=$PATH"; echo "  ALNtoPAF -> $(command -v ALNtoPAF || echo none)"' )

case_run() {  # <name> <N> <PATH for ranks> [extra env assignments...]
  local name=$1 N=$2 P=$3; shift 3; local d=$L/$name
  mkdir -p $d; cd $d
  local t0=$(date +%s.%N)
  ( export PATH=$P "$@"; $SRUN -t 3 -n $N -c 16 --cpu-bind=cores sh -c 'echo "rank $SLURM_PROCID ALNtoPAF -> $(command -v ALNtoPAF || echo none)" >&2; exec "$0" "$@"' $BIN $SMALL > merge.out 2> merge.err )
  local rc=$?
  echo "  rc=$rc  elapsed=$(echo "$(date +%s.%N) - $t0" | bc | cut -c1-5) s"
  echo "  files: $(ls | tr '\n' ' ')"
  grep -hE 'ALNtoPAF ->|ALNtoPAF failed|FAILED|\[merge\] done|limit' merge.err | sed 's/^/  log: /'
  for f in out.rank*.alntopaf.err; do [ -e "$f" ] && echo "  $f: $(head -c 200 $f | tr '\n' ' ')"; done
  cd - >/dev/null
}

echo; echo "### 2. ALNtoPAF not on PATH, N=2 (expect rc=4) ###"
case_run missing 2 /usr/bin:/bin

echo; echo "### 3. hanging fake ALNtoPAF, limit 5 s, N=2 (expect rc=4 in ~5-15 s) ###"
case_run hang 2 $FAKE:/usr/bin:/bin PGA_ALNTOPAF_TIMEOUT=5
echo "  leftover fake sleeps on node: $(pgrep -u $USER -fc 'sleep 600' || true)"
echo "### DONE $(date -Is)"
