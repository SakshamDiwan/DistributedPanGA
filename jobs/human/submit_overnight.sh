#!/bin/bash
# Overnight sweep for the report on the match-fastga-order build: the 12
# configurations of the report's timing table, three repeats each (36 jobs),
# every run compared with the ORIGINAL FastGA gold. Results are filed under
# ~/pga-report-results-order in the layout report/fig/make_figs.py reads.
set -e
cd ~/pga-jobs/human
F=$SCRATCH/pga-fixtures
export PGA_RESULTS=$HOME/pga-report-results-order; mkdir -p $PGA_RESULTS
run() { # dataset-name nodes ranks time dataset-dir stem fastga-threads
  for rep in 1 2 3; do
    sbatch --job-name=$1-n$2-r$3-rep$rep --nodes=$2 --ntasks=$3 --time=$4 \
      --export=ALL,PGA_D=$F/$5,PGA_STEM=$6,PGA_FASTGA_THREADS=$7,PGA_DS=$1,PGA_REP=$rep 07_report_sweep.sbatch
  done
}
for n in 1 2 4 8; do run yeast7     1 $n 00:15:00 yeast yeast7 8; done
for n in 4 8 16;  do run primates4  1 $n 00:40:00 primates4-chr6 primates4_chr6 32; done
run primates14 1 16 01:00:00 primates14-chr6 primates14_chr6 32
run primates14 2 32 01:00:00 primates14-chr6 primates14_chr6 32
run primates14 4 64 01:00:00 primates14-chr6 primates14_chr6 32
run hsapiens90 4 32 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32
run hsapiens90 4 64 00:40:00 hsapiens90-chr6 hsapiens90_chr6 32
