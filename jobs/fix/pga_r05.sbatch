#!/bin/bash
# R-05 check: ALNtoPAF time limit / stderr / exit code in pga-mpi-merge (MPI path).
#   1. small-fixture e2e (verify_e2e_small.sh) at N=1,2,4,8 -> must PASS, exit 0
#   2. ALNtoPAF missing from PATH, N=2        -> exit 4, .paf.partial, "cannot exec" in .err
#   3. fake hanging ALNtoPAF, limit 5 s, N=2  -> exit 4 within ~20 s, "timed out"
#   4. yeast7 N=8 with the real ALNtoPAF      -> exit 0, PAF matches gold
# Branch fix/alntopaf-timeout. Submit: sbatch ~/pga-jobs/fix/pga_r05.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=32
#SBATCH --time=00:20:00
#SBATCH --job-name=pga-r05
#SBATCH --output=%x.%j.log
set -uo pipefail
module load cpu >/dev/null 2>&1; module unload darshan >/dev/null 2>&1
R=$HOME/DistributedPanGA; cd "$R"
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1 MAKEFLAGS=CC=cc
export PATH="$HOME/FASTGA:$PATH"
L=$SCRATCH/pga-r05-$SLURM_JOB_ID; mkdir -p "$L"
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)  LOGS: $L"
echo "COMMIT: $(git rev-parse HEAD)  BRANCH: $(git rev-parse --abbrev-ref HEAD)  DIRTY: $(git status --porcelain | tr '\n' ' ')"
make -B CC=cc all > $L/build.log 2>&1; echo "build rc=$?  darshan-linked: $(ldd pga-mpi-merge | grep -c darshan)"
BIN=$L/pga-mpi-merge; cp pga-mpi-merge $BIN; sha256sum $BIN
SMALL=$R/tests/fixtures/generated/small6

echo; echo "### 1. small-fixture e2e, N=1,2,4,8 ###"
RANK_COUNTS="1 2 4 8" bash scripts/verify_e2e_small.sh > $L/e2e.log 2>&1; echo "e2e rc=$?"
grep -E 'exit=|PASS|FAIL' $L/e2e.log

case_run() {  # <name> <N> <srun --export spec> ; prints rc, elapsed, files, messages
  local name=$1 N=$2 exp=$3 d=$L/$1
  mkdir -p $d; cd $d
  local t0=$(date +%s.%N)
  srun -t 3 -n $N -c 16 --cpu-bind=cores --export=$exp $BIN $SMALL > merge.out 2> merge.err
  local rc=$?
  echo "  rc=$rc  elapsed=$(echo "$(date +%s.%N) - $t0" | bc | cut -c1-5) s"
  echo "  files: $(ls | tr '\n' ' ')"
  grep -hE 'ALNtoPAF failed|FAILED|\[merge\] done' merge.err | sed 's/^/  log: /'
  for f in out.rank*.alntopaf.err; do [ -e "$f" ] && echo "  $f: $(head -c 200 $f | tr '\n' ' ')"; done
  cd - >/dev/null
}

echo; echo "### 2. ALNtoPAF not on PATH, N=2 (expect rc=4) ###"
case_run missing 2 "ALL,PATH=/usr/bin:/bin"

echo; echo "### 3. hanging fake ALNtoPAF, limit 5 s, N=2 (expect rc=4 in ~5-15 s) ###"
FAKE=$L/fakebin; mkdir -p $FAKE
printf '#!/bin/sh\necho "fake ALNtoPAF: hanging on purpose" >&2\nsleep 600\n' > $FAKE/ALNtoPAF; chmod +x $FAKE/ALNtoPAF
case_run hang 2 "ALL,PATH=$FAKE:/usr/bin:/bin,PGA_ALNTOPAF_TIMEOUT=5"
echo "  leftover fake sleeps on node: $(pgrep -u $USER -fc 'sleep 600' || true)"

echo; echo "### 4. yeast7 N=8, real ALNtoPAF (expect rc=0, gold MATCH) ###"
Y=$R/tests/fixtures/yeast; d=$L/yeast-N8; mkdir -p $d; cd $d
t0=$(date +%s.%N)
srun -t 10 -n 8 -c 16 --cpu-bind=cores $BIN $Y/yeast7 > merge.out 2> merge.err; rc=$?
echo "  rc=$rc  wall=$(echo "$(date +%s.%N) - $t0" | bc | cut -c1-6) s  alntopaf.err files: $(ls out.rank*.alntopaf.err 2>/dev/null | wc -l)"
cat out.rank*.paf | LC_ALL=C sort | cmp -s - <(LC_ALL=C sort $Y/gold.paf) && echo "  gold=MATCH ($(cat out.rank*.paf | wc -l) alignments)" || echo "  gold=DIFFER"
cd "$R"
echo "### DONE $(date -Is)"
