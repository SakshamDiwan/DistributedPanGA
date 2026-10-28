#!/bin/bash
# Per-rank wrapper: records /usr/bin/time -v (MaxRSS, elapsed) per rank in the CWD.
exec /usr/bin/time -v -o "time.rank${SLURM_PROCID}" "$@"
