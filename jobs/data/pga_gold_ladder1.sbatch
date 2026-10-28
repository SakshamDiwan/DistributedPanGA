#!/bin/bash
# FastGA PAF gold for dataset ladder step 1-2 (Zenodo 7937947):
#   scerevisiae8 (96 Mbp), primates4.chr6 (679 Mbp), primates14.chr6 (~2.7 Gbp)
# 1. scerevisiae8 at -T8 and -T32: is FastGA's PAF thread-count invariant?
# 2. primates at -T32 if (1) is identical, else -T8.
# No pipeline binary is run. Submit: sbatch ~/pga-jobs/data/pga_gold_ladder1.sbatch
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --job-name=pga-gold-l1
#SBATCH --output=%x.%j.log
set -uo pipefail
module load cpu >/dev/null 2>&1
F=$SCRATCH/pga-fixtures; G=$HOME/pga-jobs/data/make_gold_paf.sh
echo "JOB: $SLURM_JOB_ID  NODE: $(hostname)  DATE: $(date -Is)"
sha256sum $HOME/FASTGA/{FAtoGDB,GIXmake,FastGA,ALNtoPAF}
check() { ( cd "$1" && echo "$2  $3" | md5sum -c - ) || { echo "md5 FAILED for $1/$3"; return 1; }; }
check $F/scerevisiae8   48e911dab1064fd8ca22d4e2da3d4135 scerevisiae8.fa.gz   || exit 1
check $F/primates4-chr6 60c1936d3ff940b7029b2a6ae35d0e5a primates4.chr6.fa.gz || exit 1
P14=1; check $F/primates14-chr6 5872fe857811f3c1249bccfd2b21165c primates14.chr6.fa.gz || P14=0

echo; echo "### 1. scerevisiae8: -T8 vs -T32 ###"
bash $G $F/scerevisiae8 scerevisiae8 $F/scerevisiae8/scerevisiae8.fa.gz 8  gold
bash $G $F/scerevisiae8 scerevisiae8 $F/scerevisiae8/scerevisiae8.fa.gz 32 gold.T32
if cmp -s <(LC_ALL=C sort $F/scerevisiae8/gold.paf) <(LC_ALL=C sort $F/scerevisiae8/gold.T32.paf); then
  echo "thread-invariant: sorted PAF identical at -T8 and -T32 -> using -T32"; T=32
else
  echo "NOT thread-invariant: sorted PAFs differ -> using -T8 (same as yeast7 gold)"; T=8
  diff <(LC_ALL=C sort $F/scerevisiae8/gold.paf) <(LC_ALL=C sort $F/scerevisiae8/gold.T32.paf) | head -4
fi

echo; echo "### 2. primates4.chr6 (-T$T) ###"
bash $G $F/primates4-chr6 primates4_chr6 $F/primates4-chr6/primates4.chr6.fa.gz $T gold

echo; echo "### 3. primates14.chr6 (-T$T) ###"
if [ $P14 = 1 ]; then bash $G $F/primates14-chr6 primates14_chr6 $F/primates14-chr6/primates14.chr6.fa.gz $T gold
else echo "skipped: download not verified"; fi
du -sh $F/scerevisiae8 $F/primates4-chr6 $F/primates14-chr6
echo "### DONE $(date -Is)"
