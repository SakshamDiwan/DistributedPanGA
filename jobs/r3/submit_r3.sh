#!/bin/bash
# R3 extraction/sort experiments. Run on a login node: bash jobs/r3/submit_r3.sh [bench|gold|mpi|all]
# Needs the frozen build: ~/pga-build/r3 (clone of branch r3-extract-sort, built with
#   module load cpu; module unload darshan; make CC=cc pga-mpi-merge pga-dsort-bench
#   and the MPI test binaries: make -n CC=cc test-mpi | grep -v run_mpi.sh | bash).
# Results: ~/pga-r3-results (bench: <ds>/n<N>-r<R>/*.err + bench.tsv; gold: gold/... + gold.tsv).
set -e
cd "$(dirname "$0")"
mkdir -p ~/pga-jobs/r3/logs
F=$SCRATCH/pga-fixtures
what=${1:-all}

bench() { # dataset nodes ranks time fixture-dir stem
  sbatch --job-name=r3b-$1-n$2-r$3 --nodes=$2 --ntasks=$3 --time=$4 \
    --export=ALL,PGA_D=$F/$5,PGA_STEM=$6,PGA_DS=$1 r3_dsort_bench.sbatch
}
gold() { # dataset nodes ranks time fixture-dir stem fastga-threads variant
  sbatch --job-name=r3g-$1-n$2-r$3-$8 --nodes=$2 --ntasks=$3 --time=$4 \
    --export=ALL,PGA_D=$F/$5,PGA_STEM=$6,PGA_DS=$1,PGA_FASTGA_THREADS=$7,R3_VARIANT=$8 r3_gold.sbatch
}

if [ $what = bench ] || [ $what = all ]; then
  bench yeast7     1 8   00:20:00 yeast yeast7
  bench primates14 1 16  00:45:00 primates14-chr6 primates14_chr6
  bench primates14 2 32  00:45:00 primates14-chr6 primates14_chr6
  bench primates14 4 64  00:45:00 primates14-chr6 primates14_chr6
  bench hsapiens90 4 32  01:15:00 hsapiens90-chr6 hsapiens90_chr6
  bench hsapiens90 4 64  01:00:00 hsapiens90-chr6 hsapiens90_chr6
  bench hsapiens90 4 128 01:00:00 hsapiens90-chr6 hsapiens90_chr6
  bench hsapiens90 8 128 01:00:00 hsapiens90-chr6 hsapiens90_chr6
  bench hsapiens90 8 256 01:00:00 hsapiens90-chr6 hsapiens90_chr6
fi
if [ $what = gold ] || [ $what = all ]; then
  for v in default fast; do
    gold yeast7     1 1  00:20:00 yeast yeast7 8 $v
    gold yeast7     1 8  00:15:00 yeast yeast7 8 $v
    gold primates4  1 4  00:40:00 primates4-chr6 primates4_chr6 32 $v
    gold primates4  1 16 00:30:00 primates4-chr6 primates4_chr6 32 $v
    gold primates14 4 64 01:00:00 primates14-chr6 primates14_chr6 32 $v
    gold hsapiens90 4 32 00:45:00 hsapiens90-chr6 hsapiens90_chr6 32 $v
    gold hsapiens90 4 64 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32 $v
  done
  gold hsapiens90 4 32 00:45:00 hsapiens90-chr6 hsapiens90_chr6 32 keep   # same-binary memory control
fi
if [ $what = mpi ] || [ $what = all ]; then
  sbatch r3_mpi_tests.sbatch
fi
