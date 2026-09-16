#!/bin/bash
# End-to-end check: run pga-mpi-merge at several rank counts and compare its
# AGGREGATE PAF against the independent FastGA gold for the same fixture.
#
# What is and is not evidence here:
#   - FastGA producing N alignments says nothing about this pipeline. The
#     validation is that our aggregate output MATCHES that reference, and that
#     it matches at every rank count.
#   - Duplicate PAF lines are meaningful (the same pair can align more than
#     once), so the comparison sorts but never deduplicates.
#   - An individual out.rank<N>.paf may legitimately be EMPTY: a rank can own
#     contig pairs that yield no alignment. Only the aggregate must match.
#   - pga-mpi-merge returns 0 even when ALNtoPAF fails, so a zero exit is not
#     evidence of success; the aggregate is checked for non-emptiness too.
#
# REQUIRES an active Slurm allocation. Run as a plain command, never under srun.
#
# Usage:  bash scripts/verify_e2e_small.sh
# Env:    RANK_COUNTS (default "1 2 4"), FASTGA_BIN, PGA_GDB_STEM, PGA_GOLD_PAF

set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
FIXTURES=${PGA_FIXTURES:-$ROOT/tests/fixtures/generated}
GDB_STEM=${PGA_GDB_STEM:-$FIXTURES/small6}
GOLD=${PGA_GOLD_PAF:-$FIXTURES/gold.paf}
FASTGA_BIN=${FASTGA_BIN:-$HOME/FASTGA}
RANK_COUNTS=${RANK_COUNTS:-1 2 4}

# cray-mpich demands the GTL library when the GPU modules are loaded, which
# they are by default on Perlmutter. This is a CPU-only pipeline.
export MPICH_GPU_SUPPORT_ENABLED=0
export PATH="$FASTGA_BIN:$PATH"

WORK=$(mktemp -d "${TMPDIR:-/tmp}/pga-e2e-XXXXXX")
RUNROOT=$(mktemp -d "${SCRATCH:-${TMPDIR:-/tmp}}/pga-e2e-run-XXXXXX")
cleanup() { rm -rf "$WORK" "$RUNROOT"; }
trap cleanup EXIT

for f in "$GOLD" "$GDB_STEM.1gdb"; do
    [ -f "$f" ] || { echo "ERROR: $f not found. Run scripts/make_gold_small.sh first." >&2; exit 1; }
done
command -v srun     >/dev/null || { echo "ERROR: srun not found; run inside an allocation." >&2; exit 1; }
command -v ALNtoPAF >/dev/null || { echo "ERROR: ALNtoPAF not on PATH (set FASTGA_BIN)." >&2; exit 1; }

echo "[build] pga-mpi-merge"
make -C "$ROOT" pga-mpi-merge >/dev/null || { echo "ERROR: build failed" >&2; exit 1; }

sort "$GOLD" > "$WORK/gold.sorted.paf"       # sort, never -u: duplicates matter
echo "reference: $(wc -l < "$WORK/gold.sorted.paf") alignments from FastGA"

status=0
first=""
for N in $RANK_COUNTS; do
    echo "=== world_size=$N ==="
    RUN=$RUNROOT/ws$N; mkdir -p "$RUN"; cd "$RUN"
    # Outputs are written to the process CWD by relative path, so each rank
    # count gets its own directory on a shared filesystem.
    srun -n "$N" "$ROOT/pga-mpi-merge" "$GDB_STEM" > merge.log 2>&1
    rc=$?
    nfiles=$(ls out.rank*.paf 2>/dev/null | wc -l)
    cat out.rank*.paf 2>/dev/null | sort > "$WORK/got.$N.paf"
    n=$(wc -l < "$WORK/got.$N.paf")
    echo "  exit=$rc  rank files=$nfiles (empty ones are legitimate)  aggregate=$n alignments"

    if [ "$n" -eq 0 ]; then
        echo "  FAIL: aggregate PAF is empty (exit status alone is not evidence)" >&2
        status=1
    elif cmp -s "$WORK/got.$N.paf" "$WORK/gold.sorted.paf"; then
        echo "  PASS: aggregate matches FastGA gold"
    else
        echo "  FAIL: aggregate differs from FastGA gold" >&2
        diff <(head -5 "$WORK/got.$N.paf") <(head -5 "$WORK/gold.sorted.paf") | head -10 >&2
        status=1
    fi

    if [ -z "$first" ]; then first=$N; else
        cmp -s "$WORK/got.$first.paf" "$WORK/got.$N.paf" \
            && echo "  consistency: identical to ws=$first" \
            || { echo "  FAIL: differs from ws=$first" >&2; status=1; }
    fi
    cd "$ROOT"
done

echo
[ "$status" -eq 0 ] && echo "ALL PASS: aggregate PAF matches FastGA gold at N in { $RANK_COUNTS }" \
                    || echo "SOME CHECKS FAILED" >&2
exit "$status"
