#!/bin/bash
# Step 3: FastGA reference (gold) self-alignment of hsapiens90.chr6 (needs step 2).
# Same settings as the other datasets' gold: defaults (-f10, -k40), -T32 (shown thread-invariant on scerevisiae8).
# Estimate from primates4/14 gold (cost ~ haplotypes^2 x length, fits both within 1%): ~15 h (8-30 h), ~330 GB -> 48 h limit.
# Submit: sbatch ~/pga-jobs/human/03_fastga_gold.sbatch   (or --dependency=afterok:<index job>)
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --job-name=human-gold
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
echo "### FastGA -1:gold -T32 -P on Lustre"
rm -f gold.1aln gold.paf
/usr/bin/time -f "  %e s, %M kB maxRSS" FastGA -1:gold -T32 -P$TMPDIR hsapiens90_chr6 || exit 1
echo "### ALNtoPAF"
/usr/bin/time -f "  %e s, %M kB maxRSS" ALNtoPAF -T32 gold.1aln > gold.paf || exit 1
N=$(wc -l < gold.paf); echo "gold alignments: $N"
printf "# gold (%s, job %s)\ncommand: FastGA -1:gold -T32 hsapiens90_chr6; ALNtoPAF -T32 gold.1aln > gold.paf\nalignments: %s\n%s\n\n" "$(date -u +%FT%TZ)" "$SLURM_JOB_ID" "$N" "$(sha256sum gold.paf)" >> MANIFEST.txt
rm -rf "${RUN:?}/tmp"
echo "### DONE $(date -Is)"
