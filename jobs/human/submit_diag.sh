#!/bin/bash
# Submits the two diagnostic jobs for the 1-line human DIFFER.
set -e
cd ~/pga-jobs/human
sbatch 03b_gold_T16.sbatch
sbatch --nodes=4 --ntasks=32 --time=01:00:00 --job-name=human-pilot32 04_pilot.sbatch
