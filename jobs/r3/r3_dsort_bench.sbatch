#!/bin/bash
# R3: extraction + distributed sort only (pga-dsort-bench), several variants, interleaved repeats.
# No alignment, so each run takes ~1 min even on human. Every variant must give the same kseq and
# multi digests (sorted k-mer + LCP sequence, record multiset); the job prints a digest verdict.
# Placement as the report sweep: ranks/node = ranks / nodes; -c = 256 / ranks-per-node; --cpu-bind=cores.
# Submit through submit_r3.sh (sets nodes, ntasks, PGA_D, PGA_STEM, PGA_DS).
# Output: $R3_OUT/<dataset>/n<nodes>-r<ranks>/<variant>.rep<k>.err  (+ bench.tsv summary lines)
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --job-name=r3bench
#SBATCH --output=/global/u2/h/hsc53/pga-jobs/r3/logs/%x.%j.log
set -uo pipefail
module load cpu; module unload darshan
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1
unset PGA_EXCHANGE_ROUND_BYTES PGA_SORT_THREADS PGA_SORT_PARTS PGA_STAGE2_KEEP_LOCAL PGA_EXCHANGE_LARGE_COUNT PGA_BUCKET_HIST PGA_SPLITTER PGA_SORT_HEAVY

R=${R3_BUILD:-$HOME/pga-build/r3}; BENCH=$R/pga-dsort-bench
D=${PGA_D:?}; STEM=$D/${PGA_STEM:?}; DS=${PGA_DS:?}
REPS=${R3_REPS:-3}
N=$SLURM_NTASKS; NODES=$SLURM_JOB_NUM_NODES
(( N % NODES == 0 )) || { echo "ranks ($N) not divisible by nodes ($NODES)"; exit 1; }
RPN=$(( N / NODES )); (( 256 % RPN == 0 )) || { echo "256 not divisible by ranks/node ($RPN)"; exit 1; }
C=$(( 256 / RPN )); PHYS=$(( C / 2 )); (( PHYS >= 1 )) || PHYS=1
# name|env (comma-separated KEY=VALUE; PHYS = physical cores per rank)
VARIANTS=${R3_VARIANTS:-"old|PGA_STAGE2_KEEP_LOCAL=1
take|
p1024|PGA_SORT_PARTS=1024
tphys|PGA_SORT_THREADS=PHYS
p1024tphys|PGA_SORT_PARTS=1024,PGA_SORT_THREADS=PHYS
large|PGA_EXCHANGE_LARGE_COUNT=1
r256M|PGA_EXCHANGE_ROUND_BYTES=268435456
r4G|PGA_EXCHANGE_ROUND_BYTES=4294967296
opt|PGA_SPLITTER=optimal
all|PGA_SORT_PARTS=1024,PGA_SORT_THREADS=PHYS,PGA_EXCHANGE_LARGE_COUNT=1,PGA_SPLITTER=optimal"}

# Named variant sets (sbatch --export cannot carry the commas of a list).
if [ "${R3_VARSET:-}" = heavy ]; then
VARIANTS="pre|PGA_STAGE2_KEEP_LOCAL=1,PGA_SORT_PARTS=256,PGA_EXCHANGE_LARGE_COUNT=0,PGA_SORT_HEAVY=0
noheavy|PGA_SORT_HEAVY=0
default|
phys|PGA_SORT_THREADS=PHYS"
fi
O=${R3_OUT:-$HOME/pga-r3-results}/$DS/n$NODES-r$N; mkdir -p "$O"
RUN=$SCRATCH/pga-r3/$SLURM_JOB_NAME-$SLURM_JOB_ID; mkdir -p "$RUN"; cd "$RUN"
echo "JOB $SLURM_JOB_ID  NODES $SLURM_JOB_NODELIST  DATE $(date -Is)"
echo "dataset $STEM  placement: $N ranks, $NODES nodes, $RPN/node, -c $C (physical cores/rank $PHYS), --cpu-bind=cores"
echo "toolchain: $(cc --version 2>/dev/null | head -1) | $(module -t list 2>&1 | grep -E '^(PrgEnv|cray-mpich|craype/)' | tr '\n' ' ') | CFLAGS from make -n: $(make -s -C $R -Bn CC=cc pga-mpi-merge 2>/dev/null | grep -o '\-O3.*-DLCPs' | head -1)"
echo "repo: $(git -C $R rev-parse HEAD) ($(git -C $R rev-parse --abbrev-ref HEAD))  dirty: $(git -C $R status --porcelain | grep -v '^??' | tr '\n' ' ')"
echo "BENCH: sha256 $(sha256sum < $BENCH | cut -c1-16)  darshan-linked: $(ldd $BENCH | grep -c darshan)"
[ -s $STEM.1gdb ] || { echo "GDB missing: $STEM.1gdb"; exit 1; }

for rep in $(seq 1 $REPS); do
  while IFS='|' read -r name envs; do
    [ -n "$name" ] || continue
    envs=${envs//PHYS/$PHYS}
    EXP=ALL; [ -n "$envs" ] && EXP="ALL,$envs"
    err=$O/$name.rep$rep.err
    t0=$(date +%s.%N)
    srun -n $N --ntasks-per-node=$RPN -c $C --cpu-bind=cores --kill-on-bad-exit=1 --export=$EXP \
         $BENCH $STEM 1 < /dev/null > /dev/null 2> $err   # srun must not eat the variant list on stdin
    rc=$?; wall=$(echo "$(date +%s.%N) - $t0" | bc)
    b=$(grep -h '^\[bench\] rep=' $err | head -1)
    echo "$name rep$rep rc=$rc srun_wall_s=$wall env=[$envs] ${b#\[bench\] }"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $DS $NODES $N $name $rep $rc $wall "${b#\[bench\] }" >> $O/bench.tsv
  done <<< "$VARIANTS"
done

echo "--- digests (must be one value each) ---"
for k in kseq multi order; do
  echo "$k: $(grep -ho "$k=[0-9a-f]*" $O/*.err | sort | uniq -c | tr '\n' ' ')"
done
n_kseq=$(grep -ho 'kseq=[0-9a-f]*' $O/*.err | sort -u | wc -l)
n_multi=$(grep -ho 'multi=[0-9a-f]*' $O/*.err | sort -u | wc -l)
n_runs=$(grep -l '^\[bench\] rep=' $O/*.err | wc -l)
[ "$n_kseq" = 1 ] && [ "$n_multi" = 1 ] && echo "DIGESTS: SAME across $n_runs runs" || echo "DIGESTS: DIFFER (kseq $n_kseq, multi $n_multi values)"
rm -rf "${RUN:?}"
echo "### DONE $(date -Is)"
