#!/bin/bash
# Reruns 64 and 32 ranks with per-task records, to trace the missing self-alignment.
set -e
cd ~/pga-jobs/human
sbatch --nodes=4 --ntasks=64 --job-name=human-trace64 04c_trace.sbatch
sbatch --nodes=4 --ntasks=32 --job-name=human-trace32 04c_trace.sbatch
