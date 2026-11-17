#!/bin/bash
# Step 7: report sweep, match-fastga-order build (self pairs start from the value FastGA would have had) against the ORIGINAL FastGA gold.
# Same commands as 04_pilot.sbatch; only the binary differs. Prints gold=MATCH or gold=DIFFER and keeps the [fga-self] lines.
# Placement: ranks/node = ranks / nodes; -c = 256 / ranks-per-node; --cpu-bind=cores.
# Submit through submit_overnight.sh. With PGA_RESULTS, PGA_DS and PGA_REP set, also files the run under
# $PGA_RESULTS/<dataset>/n<nodes>-r<ranks>/rep<k>/ and appends a row to $PGA_RESULTS/results.tsv, the layout
# report/fig/make_figs.py reads.
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --nodes=4
#SBATCH --ntasks=64
#SBATCH --time=01:00:00
#SBATCH --job-name=order
#SBATCH --output=/global/u2/h/hsc53/pga-jobs/human/logs/%x.%j.log
set -uo pipefail
export FASTGA_DIR=$HOME/FASTGA                     # FastGA 178afcb (verified by rebuild)
export BIN=$HOME/pga-build/match-fastga/pga-mpi-merge   # branch match-fastga-order

# --- identical for both of you ---
module load cpu; module unload darshan
export PATH=$FASTGA_DIR:$PATH
export MPICH_GPU_SUPPORT_ENABLED=0
export OMP_NUM_THREADS=1
export PGA_ALNTOPAF_TIMEOUT=3600
unset  PGA_EXCHANGE_ROUND_BYTES PGA_PAIR_TASKS PGA_MEM_STAGE_RESET
export PGA_FASTGA_ORDER_LOG=1      # log every self pair entry value and the FastGA ranges
echo "PGA_FASTGA_ORDER=${PGA_FASTGA_ORDER:-1} PGA_FASTGA_THREADS=${PGA_FASTGA_THREADS:-32} PGA_FASTGA_UNSET_BEPOS=${PGA_FASTGA_UNSET_BEPOS:-32767}"

# --- run folder: fresh, on Lustre, never /tmp ---
RUN=$SCRATCH/pga-human/$SLURM_JOB_NAME-$SLURM_JOB_ID; mkdir -p "${RUN:?}/tmp"
export TMPDIR=$RUN/tmp          # sort and anything else that defaults to /tmp
# Dataset: human by default; override with --export=ALL,PGA_D=<fixture dir>,PGA_STEM=<gdb stem>
D=${PGA_D:-$SCRATCH/pga-fixtures/hsapiens90-chr6}; STEM=$D/${PGA_STEM:-hsapiens90_chr6}
echo "dataset: $STEM"
echo "JOB: $SLURM_JOB_ID  NODES: $SLURM_JOB_NODELIST  DATE: $(date -Is)  RUN: $RUN"
N=$SLURM_NTASKS; NODES=$SLURM_JOB_NUM_NODES
(( N % NODES == 0 )) || { echo "ranks ($N) not divisible by nodes ($NODES)"; exit 1; }
RPN=$(( N / NODES )); (( 256 % RPN == 0 )) || { echo "256 not divisible by ranks/node ($RPN)"; exit 1; }
C=$(( 256 / RPN ))
echo "placement: $N ranks, $NODES nodes, $RPN ranks/node, -c $C, --cpu-bind=cores"
[ -s $STEM.1gdb ] || { echo "GDB missing: run 02_index first"; exit 1; }
R=$HOME/pga-build/match-fastga
echo "repo: $(git -C $R rev-parse HEAD) ($(git -C $R rev-parse --abbrev-ref HEAD))  dirty: $(git -C $R status --porcelain | grep -v fixtures | tr '\n' ' ')"
echo "BIN: $BIN  sha256 $(sha256sum < $BIN | cut -c1-16)  built $(stat -c %y $BIN | cut -c1-19)  darshan-linked: $(ldd $BIN | grep -c darshan)"
export PGA_SCRATCH_DIR=$RUN/scratch; mkdir -p "${PGA_SCRATCH_DIR:?}"     # fresh, on Lustre
mkdir -p "${RUN:?}/out"; cd "$RUN/out"
t0=$(date +%s.%N)
srun -n $N --ntasks-per-node=$RPN -c $C --cpu-bind=cores $HOME/pga-jobs/week1/pga_w1_rankwrap.sh $BIN $STEM > merge.out 2> merge.err
rc=$?; wall=$(echo "$(date +%s.%N) - $t0" | bc)
echo "pga-mpi-merge rc=$rc wall_s=$wall"
tail -1 pga-mpi-merge.timing.csv 2>/dev/null | sed 's/^/csv: /'
grep -h '^\[pr1-mem\]' merge.err | grep -oE 'overall_peak_kb=[0-9]+' | sort -t= -k2 -n | tail -1 | sed 's/^/max rank /'
cat out.rank*.paf | LC_ALL=C sort -S 32G --parallel=16 -T $TMPDIR > got.sorted.paf
echo "alignments: $(wc -l < got.sorted.paf)"
if [ -s $D/gold.paf ]; then
  LC_ALL=C sort -S 32G --parallel=16 -T $TMPDIR $D/gold.paf > $TMPDIR/gold.sorted.paf
  cmp -s got.sorted.paf $TMPDIR/gold.sorted.paf && echo "gold=MATCH ($(wc -l < $TMPDIR/gold.sorted.paf))" || { echo "gold=DIFFER only-got $(LC_ALL=C comm -23 got.sorted.paf $TMPDIR/gold.sorted.paf | wc -l) only-gold $(LC_ALL=C comm -13 got.sorted.paf $TMPDIR/gold.sorted.paf | wc -l)"; LC_ALL=C comm -3 got.sorted.paf $TMPDIR/gold.sorted.paf | cut -c1-150 | head -6; }
else
  echo "gold not ready yet: kept $RUN/out/got.sorted.paf for a later cmp"
fi
grep -ho '\[fga-order\] on.*' merge.err | head -1
grep -ho '\[fga-self\].*' merge.err | LC_ALL=C sort -t= -k3 -n > fga-self.txt
grep -ho '\[fga-range\].*' merge.err | LC_ALL=C sort -t= -k2 -n > fga-range.txt
echo "self pairs: $(wc -l < fga-self.txt)  from earlier row: $(grep -c 'from=row' fga-self.txt)  unset: $(grep -c 'from=unset' fga-self.txt)  ranges: $(wc -l < fga-range.txt)"
grep -ho '\[fga-order\] rank.*' merge.err | awk '{for(i=1;i<=NF;i++){split($i,a,"=");if(a[1]=="wait_s"&&a[2]>w)w=a[2];if(a[1]=="rounds_s"&&a[2]>r)r=a[2];if(a[1]=="set_aside")s+=a[2];if(a[1]=="answers_sent")m+=a[2]}}END{printf "order overhead: self pairs set aside %d; longest wait for a result from another rank %.3f s; results sent between ranks %d; longest time in final rounds %.3f s\n",s,w,m,r}'
python3 -I $R/scripts/analyze_align_balance.py $RUN/out > balance.txt 2>&1; grep -E '^max/mean' balance.txt
if [ -n "${PGA_RESULTS:-}" ]; then
  DS=${PGA_DS:?}; REP=${PGA_REP:?}
  O=$PGA_RESULTS/$DS/n$NODES-r$N/rep$REP; mkdir -p "$O"
  cp merge.err merge.out pga-mpi-merge.timing.csv balance.txt fga-self.txt "$O"/ 2>/dev/null
  G=DIFFER; cmp -s got.sorted.paf $TMPDIR/gold.sorted.paf && G=MATCH
  IFS=, read -r w e s m a t < <(tail -1 pga-mpi-merge.timing.csv)
  P=$(grep -h '^\[pr1-mem\]' merge.err | grep -oE 'overall_peak_kb=[0-9]+' | cut -d= -f2 | sort -n | tail -1)
  ( flock 9
    [ -s $PGA_RESULTS/results.tsv ] || printf 'dataset\tnodes\tranks\trep\trc\textract_s\tsort_s\tmerge_s\talignments_s\ttotal_s\twall_s\tgold\tmax_rank_peak_kb\tjob\n' >&9
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $DS $NODES $N $REP $rc $e $s $m $a $t $wall $G ${P:-NA} $SLURM_JOB_ID >&9
  ) 9>>$PGA_RESULTS/results.tsv
  echo "filed: $O  ($DS n$NODES r$N rep$REP gold=$G)"
fi
rm -rf "${RUN:?}/tmp" "${RUN:?}/scratch"
sacct -j $SLURM_JOB_ID --format=JobID%-16,Elapsed,MaxRSS,State -P | tail -n +2 | head -5
echo "### DONE $(date -Is)"
