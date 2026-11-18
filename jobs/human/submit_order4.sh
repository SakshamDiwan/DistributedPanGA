#!/bin/bash
# Fourth round for branch match-fastga-order (results sent rank to rank):
# the 12 runs against the ORIGINAL FastGA gold, plus 3 control runs of the
# same binary with the feature off (PGA_FASTGA_ORDER=0) for a same-day time
# comparison. The human 64-rank control is expected to print gold=DIFFER:
# it is the old behaviour that loses one alignment.
set -e
cd ~/pga-jobs/human
F=$SCRATCH/pga-fixtures
run() { # name nodes ranks time dataset-dir stem fastga-threads order(1/0)
  sbatch --job-name=$1 --nodes=$2 --ntasks=$3 --time=$4 \
    --export=ALL,PGA_D=$F/$5,PGA_STEM=$6,PGA_FASTGA_THREADS=$7,PGA_FASTGA_ORDER=${8:-1} 06_order.sbatch
}
for n in 1 2 4 8;  do run yeast7-or4-$n 1 $n 00:15:00 yeast yeast7 8; done
for n in 4 8 16;   do run prim4-or4-$n  1 $n 00:40:00 primates4-chr6 primates4_chr6 32; done
run prim14-or4-16 1 16 01:00:00 primates14-chr6 primates14_chr6 32
run prim14-or4-32 2 32 01:00:00 primates14-chr6 primates14_chr6 32
run prim14-or4-64 4 64 01:00:00 primates14-chr6 primates14_chr6 32
run human-or4-32  4 32 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32
run human-or4-64  4 64 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32
run prim14-off-64 4 64 01:00:00 primates14-chr6 primates14_chr6 32 0
run human-off-32  4 32 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32 0
run human-off-64  4 64 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32 0
