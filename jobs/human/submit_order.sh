#!/bin/bash
# Tests branch match-fastga-order against the ORIGINAL FastGA gold at several
# rank counts per dataset (12 runs), plus the diagnostic FastGA on three
# datasets (3 runs) to check the entry values directly.
set -e
cd ~/pga-jobs/human
F=$SCRATCH/pga-fixtures
run() { # name nodes ranks time dataset-dir stem fastga-threads
  sbatch --job-name=$1 --nodes=$2 --ntasks=$3 --time=$4 \
    --export=ALL,PGA_D=$F/$5,PGA_STEM=$6,PGA_FASTGA_THREADS=$7 06_order.sbatch
}
for n in 1 2 4 8;  do run yeast7-ord$n 1 $n 00:15:00 yeast yeast7 8; done
for n in 4 8 16;   do run prim4-ord$n  1 $n 00:40:00 primates4-chr6 primates4_chr6 32; done
run prim14-ord16 1 16 01:30:00 primates14-chr6 primates14_chr6 32
run prim14-ord32 2 32 01:00:00 primates14-chr6 primates14_chr6 32
run prim14-ord64 4 64 01:00:00 primates14-chr6 primates14_chr6 32
run human-ord32  4 32 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32
run human-ord64  4 64 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32
diag() { sbatch --job-name=$1 --time=$2 --export=ALL,PGA_D=$F/$3,PGA_STEM=$4,PGA_T=32 06b_fastga_diag.sbatch; }
diag prim4-diag  00:20:00 primates4-chr6 primates4_chr6
diag prim14-diag 01:00:00 primates14-chr6 primates14_chr6
diag human-diag  01:00:00 hsapiens90-chr6 hsapiens90_chr6
