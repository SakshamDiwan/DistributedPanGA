#!/bin/bash
# P5 runner. Must be run INSIDE an allocation, as a plain command -- never
# under srun, since it calls srun itself.
set -uo pipefail

# The default Perlmutter login/compute environment loads craype-accel-nvidia80
# and gpu/1.0, so cray-mpich demands the GTL library and aborts with
# "GPU_SUPPORT_ENABLED is requested, but GTL library is not linked".
# This is a CPU-only project, so turn that off. (`module load cpu` also works.)
export MPICH_GPU_SUPPORT_ENABLED=0

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=$ROOT/tests/bin/test_stage2
WORK=${TMPDIR:-/tmp}/pga-mpi-tests.$$
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

if ! command -v srun >/dev/null 2>&1; then
    echo "ERROR: srun not found. Run inside a Slurm allocation." >&2
    exit 1
fi
if [ ! -x "$BIN" ]; then
    echo "ERROR: $BIN not built. Run 'make test-mpi'." >&2
    exit 1
fi

RANK_COUNTS="${RANK_COUNTS:-1 2 3 4 8}"
status=0

for N in $RANK_COUNTS; do
    echo "=== world_size=$N ==="
    srun -n "$N" "$BIN" uniform "$WORK/digest.$N" || status=1
    srun -n "$N" "$BIN" skew                      || status=1
    srun -n "$N" "$BIN" allempty                  || status=1
done

# R-01 needs exactly 4 ranks.
case " $RANK_COUNTS " in
    *" 4 "*) echo "=== R-01 seam (world_size=4) ==="
             srun -n 4 "$BIN" seam || status=1 ;;
    *)       echo "note: R-01 seam skipped (RANK_COUNTS has no 4)" ;;
esac

# M-06: the canonical digest must be identical at every world size.
echo "=== M-06 cross-world-size consistency ==="
ref=""
for N in $RANK_COUNTS; do
    [ -f "$WORK/digest.$N" ] || continue
    if [ -z "$ref" ]; then ref=$WORK/digest.$N; refn=$N; continue; fi
    if cmp -s "$ref" "$WORK/digest.$N"; then
        echo "  ws=$N matches ws=$refn"
    else
        echo "  FAIL: ws=$N differs from ws=$refn" >&2
        status=1
    fi
done

[ "$status" -eq 0 ] && echo "ALL MPI TESTS PASS" || echo "SOME MPI TESTS FAILED" >&2
exit "$status"
