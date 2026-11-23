#!/bin/bash
# R3: end-to-end gold check of the r3-extract-sort build (pga-mpi-merge) against the ORIGINAL FastGA gold.
# Same commands as jobs/human/07_report_sweep.sbatch; only the binary and the optional sort/exchange
# environment differ. R3_VARIANT names the environment (default = build defaults, which now free the
# extraction buffer after packing). Prints gold=MATCH or gold=DIFFER.
# Files results under $R3_OUT/gold/<dataset>/n<nodes>-r<ranks>/<variant>/ and appends to gold.tsv.
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --job-name=r3gold
#SBATCH --output=/global/u2/h/hsc53/pga-jobs/r3/logs/%x.%j.log
set -uo pipefail
export FASTGA_DIR=$HOME/FASTGA                      # FastGA 178afcb
R=$HOME/pga-build/r3; export BIN=$R/pga-mpi-merge
module load cpu; module unload darshan
export PATH=$FASTGA_DIR:$PATH
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1
export PGA_ALNTOPAF_TIMEOUT=3600
unset  PGA_EXCHANGE_ROUND_BYTES PGA_PAIR_TASKS PGA_MEM_STAGE_RESET PGA_SORT_THREADS PGA_SORT_PARTS PGA_STAGE2_KEEP_LOCAL PGA_EXCHANGE_LARGE_COUNT

N=$SLURM_NTASKS; NODES=$SLURM_JOB_NUM_NODES
(( N % NODES == 0 )) || { echo "ranks ($N) not divisible by nodes ($NODES)"; exit 1; }
RPN=$(( N / NODES )); (( 256 % RPN == 0 )) || { echo "256 not divisible by ranks/node ($RPN)"; exit 1; }
C=$(( 256 / RPN )); PHYS=$(( C / 2 )); (( PHYS >= 1 )) || PHYS=1
VARIANT=${R3_VARIANT:-default}
case $VARIANT in
  default) ;;
  keep)    export PGA_STAGE2_KEEP_LOCAL=1 ;;
  fast)    export PGA_SORT_PARTS=1024 PGA_SORT_THREADS=$PHYS PGA_EXCHANGE_LARGE_COUNT=1 ;;
  *) echo "unknown R3_VARIANT $VARIANT"; exit 1 ;;
esac
echo "variant $VARIANT: PGA_SORT_PARTS=${PGA_SORT_PARTS:-} PGA_SORT_THREADS=${PGA_SORT_THREADS:-} PGA_EXCHANGE_LARGE_COUNT=${PGA_EXCHANGE_LARGE_COUNT:-} PGA_STAGE2_KEEP_LOCAL=${PGA_STAGE2_KEEP_LOCAL:-}"
echo "PGA_FASTGA_ORDER=${PGA_FASTGA_ORDER:-1} PGA_FASTGA_THREADS=${PGA_FASTGA_THREADS:-32}"

RUN=$SCRATCH/pga-r3/$SLURM_JOB_NAME-$SLURM_JOB_ID; mkdir -p "${RUN:?}/tmp"
export TMPDIR=$RUN/tmp
D=${PGA_D:?}; STEM=$D/${PGA_STEM:?}; DS=${PGA_DS:?}
echo "dataset: $STEM"
echo "JOB: $SLURM_JOB_ID  NODES: $SLURM_JOB_NODELIST  DATE: $(date -Is)  RUN: $RUN"
echo "placement: $N ranks, $NODES nodes, $RPN ranks/node, -c $C, --cpu-bind=cores"
[ -s $STEM.1gdb ] || { echo "GDB missing"; exit 1; }
[ -s $D/gold.paf ] || { echo "gold missing: $D/gold.paf"; exit 1; }
echo "toolchain: $(cc --version 2>/dev/null | head -1) | $(module -t list 2>&1 | grep -E '^(PrgEnv|cray-mpich|craype/)' | tr '\n' ' ') | CFLAGS from make -n: $(make -s -C $R -Bn CC=cc pga-mpi-merge 2>/dev/null | grep -o '\-O3.*-DLCPs' | head -1)"
echo "repo: $(git -C $R rev-parse HEAD) ($(git -C $R rev-parse --abbrev-ref HEAD))  dirty: $(git -C $R status --porcelain | grep -v '^??' | tr '\n' ' ')"
echo "BIN: $BIN  sha256 $(sha256sum < $BIN | cut -c1-16)  built $(stat -c %y $BIN | cut -c1-19)  darshan-linked: $(ldd $BIN | grep -c darshan)"
export PGA_SCRATCH_DIR=$RUN/scratch; mkdir -p "${PGA_SCRATCH_DIR:?}"
mkdir -p "${RUN:?}/out"; cd "$RUN/out"
t0=$(date +%s.%N)
srun -n $N --ntasks-per-node=$RPN -c $C --cpu-bind=cores $HOME/pga-jobs/week1/pga_w1_rankwrap.sh $BIN $STEM > merge.out 2> merge.err
rc=$?; wall=$(echo "$(date +%s.%N) - $t0" | bc)
echo "pga-mpi-merge rc=$rc wall_s=$wall"
tail -1 pga-mpi-merge.timing.csv 2>/dev/null | sed 's/^/csv: /'
P=$(grep -h '^\[pr1-mem\]' merge.err | grep -oE 'overall_peak_kb=[0-9]+' | cut -d= -f2 | sort -n | tail -1)
echo "max rank overall_peak_kb=$P"
echo "[r3-s2] lines: $(grep -c '^\[r3-s2\]' merge.err)  took_local: $(grep -ho 'took_local=[01]' merge.err | sort | uniq -c | tr '\n' ' ')  sort: $(grep -ho 'sort_threads=[0-9]* sort_parts=[0-9]*' merge.err | sort | uniq -c | tr '\n' ' ')"
cat out.rank*.paf | LC_ALL=C sort -S 32G --parallel=16 -T $TMPDIR > got.sorted.paf
LC_ALL=C sort -S 32G --parallel=16 -T $TMPDIR $D/gold.paf > $TMPDIR/gold.sorted.paf
echo "alignments: $(wc -l < got.sorted.paf)  gold: $(wc -l < $TMPDIR/gold.sorted.paf)"
G=DIFFER
if cmp -s got.sorted.paf $TMPDIR/gold.sorted.paf; then G=MATCH; echo "gold=MATCH ($(wc -l < $TMPDIR/gold.sorted.paf))"
else
  echo "gold=DIFFER only-got $(LC_ALL=C comm -23 got.sorted.paf $TMPDIR/gold.sorted.paf | wc -l) only-gold $(LC_ALL=C comm -13 got.sorted.paf $TMPDIR/gold.sorted.paf | wc -l)"
fi
O=${R3_OUT:-$HOME/pga-r3-results}/gold/$DS/n$NODES-r$N/$VARIANT; mkdir -p "$O"
cp merge.err merge.out pga-mpi-merge.timing.csv "$O"/ 2>/dev/null
IFS=, read -r w e s m a t < <(tail -1 pga-mpi-merge.timing.csv)
( flock 9
  [ -s ${R3_OUT:-$HOME/pga-r3-results}/gold.tsv ] || printf 'dataset\tnodes\tranks\tvariant\trc\textract_s\tsort_s\tmerge_s\talignments_s\ttotal_s\twall_s\tgold\tmax_rank_peak_kb\tjob\n' >&9
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $DS $NODES $N $VARIANT $rc $e $s $m $a $t $wall $G ${P:-NA} $SLURM_JOB_ID >&9
) 9>>${R3_OUT:-$HOME/pga-r3-results}/gold.tsv
echo "filed: $O ($DS n$NODES r$N $VARIANT gold=$G)"
rm -rf "${RUN:?}"
sacct -j $SLURM_JOB_ID --format=JobID%-16,Elapsed,MaxRSS,State -P | tail -n +2 | head -5
echo "### DONE $(date -Is)"
