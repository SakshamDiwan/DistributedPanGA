#!/bin/bash
# Second round for branch match-fastga-order (reworked schedule): the same 12
# runs against the ORIGINAL FastGA gold. No diagnostic FastGA runs this time.
set -e
cd ~/pga-jobs/human
F=$SCRATCH/pga-fixtures
run() { # name nodes ranks time dataset-dir stem fastga-threads
  sbatch --job-name=$1 --nodes=$2 --ntasks=$3 --time=$4 \
    --export=ALL,PGA_D=$F/$5,PGA_STEM=$6,PGA_FASTGA_THREADS=$7 06_order.sbatch
}
for n in 1 2 4 8;  do run yeast7-or2-$n 1 $n 00:15:00 yeast yeast7 8; done
for n in 4 8 16;   do run prim4-or2-$n  1 $n 00:40:00 primates4-chr6 primates4_chr6 32; done
run prim14-or2-16 1 16 01:00:00 primates14-chr6 primates14_chr6 32
run prim14-or2-32 2 32 01:00:00 primates14-chr6 primates14_chr6 32
run prim14-or2-64 4 64 01:00:00 primates14-chr6 primates14_chr6 32
run human-or2-32  4 32 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32
run human-or2-64  4 64 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32
