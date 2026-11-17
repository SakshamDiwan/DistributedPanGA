#!/bin/bash
# Third round for branch match-fastga-order (unset self pairs resolved at once):
# timing check on the two large datasets, plus one small run of each other
# dataset as a regression check. All against the ORIGINAL FastGA gold.
set -e
cd ~/pga-jobs/human
F=$SCRATCH/pga-fixtures
run() { # name nodes ranks time dataset-dir stem fastga-threads
  sbatch --job-name=$1 --nodes=$2 --ntasks=$3 --time=$4 \
    --export=ALL,PGA_D=$F/$5,PGA_STEM=$6,PGA_FASTGA_THREADS=$7 06_order.sbatch
}
run yeast7-or3-8  1 8  00:15:00 yeast yeast7 8
run prim4-or3-16  1 16 00:40:00 primates4-chr6 primates4_chr6 32
run prim14-or3-16 1 16 01:00:00 primates14-chr6 primates14_chr6 32
run prim14-or3-32 2 32 01:00:00 primates14-chr6 primates14_chr6 32
run prim14-or3-64 4 64 01:00:00 primates14-chr6 primates14_chr6 32
run human-or3-32  4 32 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32
run human-or3-64  4 64 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32
