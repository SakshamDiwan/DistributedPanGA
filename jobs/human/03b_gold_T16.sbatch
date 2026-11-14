#!/bin/bash
# Step 3b: FastGA rerun at -T16 to test thread-invariance on human (gold -T32 had 1 line pilot lacked). Writes to the run dir; leaves gold.paf alone.
# Same settings as the other datasets' gold: defaults (-f10, -k40), -T32 (shown thread-invariant on scerevisiae8).
# Submit: sbatch ~/pga-jobs/human/03_fastga_gold.sbatch   (or --dependency=afterok:<index job>)
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --job-name=human-gold16
#SBATCH --output=%x.%j.log
set -uo pipefail
export FASTGA_DIR=$HOME/FASTGA                     # FastGA 178afcb (verified by rebuild)
export BIN=$HOME/DistributedPanGA/pga-mpi-merge    # your own build (same commit, same command)

# --- identical for both of you ---
module load cpu; module unload darshan
export PATH=$FASTGA_DIR:$PATH
export MPICH_GPU_SUPPORT_ENABLED=0
export OMP_NUM_THREADS=1
export PGA_ALNTOPAF_TIMEOUT=3600
unset  PGA_EXCHANGE_ROUND_BYTES PGA_PAIR_TASKS PGA_MEM_STAGE_RESET

# --- run folder: fresh, on Lustre, never /tmp ---
RUN=$SCRATCH/pga-human/$SLURM_JOB_NAME-$SLURM_JOB_ID; mkdir -p "${RUN:?}/tmp"
export TMPDIR=$RUN/tmp          # sort and anything else that defaults to /tmp
D=$SCRATCH/pga-fixtures/hsapiens90-chr6; STEM=$D/hsapiens90_chr6
echo "JOB: $SLURM_JOB_ID  NODES: $SLURM_JOB_NODELIST  DATE: $(date -Is)  RUN: $RUN"
cd $D
[ -s hsapiens90_chr6.gix ] || { echo "index missing: run 02_index first"; exit 1; }
echo "### FastGA -T16 -> $RUN/gold16"
/usr/bin/time -f "  %e s, %M kB maxRSS" FastGA -1:$RUN/gold16 -T16 -P$TMPDIR hsapiens90_chr6 || exit 1
echo "### ALNtoPAF"
/usr/bin/time -f "  %e s, %M kB maxRSS" ALNtoPAF -T16 $RUN/gold16.1aln > $RUN/gold16.paf || exit 1
N=$(wc -l < $RUN/gold16.paf); echo "gold16 alignments: $N (gold -T32: $(wc -l < gold.paf))"
LC_ALL=C sort -S 32G --parallel=16 -T $TMPDIR $RUN/gold16.paf > $RUN/gold16.sorted.paf
LC_ALL=C sort -S 32G --parallel=16 -T $TMPDIR gold.paf > $TMPDIR/gold32.sorted.paf
cmp -s $RUN/gold16.sorted.paf $TMPDIR/gold32.sorted.paf && echo "T16 vs T32: MATCH" || { echo "T16 vs T32: DIFFER"; LC_ALL=C comm -3 $RUN/gold16.sorted.paf $TMPDIR/gold32.sorted.paf | head -20; }
cmp -s $RUN/gold16.sorted.paf $SCRATCH/pga-human/human-pilot-61244457/out/got.sorted.paf && echo "T16 vs pilot64: MATCH" || echo "T16 vs pilot64: DIFFER"
rm -rf "${RUN:?}/tmp"
echo "### DONE $(date -Is)"
