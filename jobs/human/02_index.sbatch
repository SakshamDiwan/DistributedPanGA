#!/bin/bash
# Step 2: build the FastGA GDB + index for hsapiens90.chr6 (needs step 1).
# Submit: sbatch ~/pga-jobs/human/02_index.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --job-name=human-index
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
echo "a6513914ade700e4aef4a29290f4bec8  hsapiens90.chr6.masked.fa.gz" | md5sum -c - || exit 1
[ -e hsapiens90_chr6.fa.gz ] || ln -s hsapiens90.chr6.masked.fa.gz hsapiens90_chr6.fa.gz
sha256sum $FASTGA_DIR/{FAtoGDB,GIXmake,FastGA,ALNtoPAF} | sed 's|/.*/||'
echo "### input summary"
zcat hsapiens90_chr6.fa.gz | awk '/^>/{n++; next} {L=length($0); t+=L; lc+=gsub(/[acgtn]/,""); N+=gsub(/[Nn]/,"")}
  END{printf "  sequences=%d  bases=%.2f Gbp  lowercase(soft-masked)=%.2f%%  N=%.2f%%\n", n, t/1e9, 100*lc/t, 100*N/t}'
# Keep an existing GDB: the pilot may be reading it while this job runs.
if [ -s hsapiens90_chr6.1gdb ] && [ -s .hsapiens90_chr6.bps ]; then
  echo "### FAtoGDB: skipped, GDB already built"; ls -la hsapiens90_chr6.1gdb .hsapiens90_chr6.bps
else
  echo "### FAtoGDB"; rm -f hsapiens90_chr6.1gdb
  /usr/bin/time -f "  %e s, %M kB maxRSS" FAtoGDB hsapiens90_chr6.fa.gz || exit 1
fi
# -T16, not -T32: FastGA 178afcb GIXmake.c:2001 computes its work-block size in
# 32-bit int, ((SCAN_MAX+8) + NPARTS*BUFF_MAX)*NTHREADS. Human chr6 gives NPARTS=64,
# so at -T32 that is 2.37e9 -> overflows, the block is under-allocated and push()
# segfaults (job 61239144, signal 11 after 6 s). Any -T <= 29 fits; 16 is safe.
echo "### GIXmake -T16 -P on Lustre"
/usr/bin/time -f "  %e s, %M kB maxRSS" GIXmake -T16 -P$TMPDIR hsapiens90_chr6 || exit 1
ls -la hsapiens90_chr6.1gdb hsapiens90_chr6.gix; du -sh $D
printf "# index (%s, job %s)\nfastga:  178afcb at %s\ncommand: FAtoGDB hsapiens90_chr6.fa.gz; GIXmake -T16 -P<run tmp> hsapiens90_chr6\n\n" "$(date -u +%FT%TZ)" "$SLURM_JOB_ID" "$FASTGA_DIR" >> MANIFEST.txt
rm -rf "${RUN:?}/tmp"
echo "### DONE $(date -Is)"
