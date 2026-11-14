#!/bin/bash
# Logs every aligner call for contig 1072 at 64 and 32 ranks (debug build).
set -e
cd ~/pga-jobs/human
sbatch --nodes=4 --ntasks=64 --job-name=human-ta64 04d_trace_align.sbatch
sbatch --nodes=4 --ntasks=32 --job-name=human-ta32 04d_trace_align.sbatch
