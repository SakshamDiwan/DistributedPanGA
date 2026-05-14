#!/bin/bash
# Phase 4b runtime test driver.
# See src/phase4b_runtime_testing_plan.md.
#
# Behavior depends on $SLURM_JOB_NUM_NODES of the enclosing allocation:
#   1 node:  single-node sweep at world_size 1, 2, 4, 8 + correctness gate
#            + final stage × world_size table from pga-mpi-merge.timing.csv.
#   2 nodes: cross-node correctness at 2x1 (total=2) and 2x4 (total=8).
#   4 nodes: cross-node correctness at 4x2 (total=8) and 4x4 (total=16).
#
# Usage:
#   salloc -N 1 -C cpu -q debug -t 30
#   bash scripts/runtime_test_phase4b.sh
#
# Logs land in $ROOT/logs/. CSV at $ROOT/pga-mpi-merge.timing.csv (appended).

set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
GDB_STEM=$ROOT/tests/fixtures/yeast/yeast7
GOLD=$ROOT/tests/fixtures/yeast/gold.paf
FASTGA_BIN_DIR=/global/u2/s/sdiwan/FASTGA

if [ ! -f "$GOLD" ]; then
    echo "ERROR: $GOLD not found. Run scripts/make_gold.sh first." >&2
    exit 1
fi
if [ ! -f "$GDB_STEM.1gdb" ]; then
    echo "ERROR: $GDB_STEM.1gdb not found." >&2
    exit 1
fi
if ! command -v srun >/dev/null 2>&1; then
    echo "ERROR: srun not found. Run inside a Slurm allocation." >&2
    exit 1
fi
if [ ! -x "$FASTGA_BIN_DIR/ALNtoPAF" ]; then
    echo "ERROR: $FASTGA_BIN_DIR/ALNtoPAF not found." >&2
    exit 1
fi

echo "[build] pga-mpi-merge..."
make -C "$ROOT" pga-mpi-merge >/dev/null

mkdir -p "$ROOT/logs"
sort "$GOLD" > /tmp/gold.sorted.paf

# Truncate the CSV at the start of every run so the table reflects only this
# invocation's measurements.
rm -f "$ROOT/pga-mpi-merge.timing.csv"

run_one() {
    # $1=label, rest=srun args before the binary path
    local label="$1"; shift
    rm -f "$ROOT"/out.rank*.paf
    PATH=$FASTGA_BIN_DIR:$PATH \
        srun "$@" "$ROOT/pga-mpi-merge" "$GDB_STEM" \
        2>&1 | tee "$ROOT/logs/timing.${label}.log"
    cat "$ROOT"/out.rank*.paf | sort > /tmp/got.sorted.paf
    if cmp -s /tmp/got.sorted.paf /tmp/gold.sorted.paf; then
        echo "  $label CORRECTNESS PASS"
        return 0
    else
        echo "  $label CORRECTNESS FAIL — see $ROOT/logs/timing.${label}.log" >&2
        return 1
    fi
}

NODES=${SLURM_JOB_NUM_NODES:-1}

case "$NODES" in
    1)
        for N in 1 2 4 8; do
            echo "=== single-node N=$N ==="
            run_one "singlenode.N$N" -n "$N"
        done

        # After the sweep, render the stage × world_size table.
        if [ -s "$ROOT/pga-mpi-merge.timing.csv" ]; then
            echo
            echo "=== Stage × world_size (max across ranks, seconds) ==="
            awk -F, '
                NR == 1 { next }                      # header
                {
                    ws = $1
                    extract_s[ws] = $2
                    sort_s[ws]    = $3
                    merge_s[ws]   = $4
                    align_s[ws]   = $5
                    total_s[ws]   = $6
                    seen[ws] = 1
                }
                END {
                    n = asorti(seen, ws_sorted, "@val_num_asc")
                    printf "%-12s", "Stage"
                    for (i = 1; i <= n; i++) printf "  ws=%-7s", ws_sorted[i]
                    printf "\n"

                    rows[1] = "extract";    vals[1] = "extract_s"
                    rows[2] = "sort";       vals[2] = "sort_s"
                    rows[3] = "merge";      vals[3] = "merge_s"
                    rows[4] = "alignments"; vals[4] = "align_s"
                    rows[5] = "TOTAL";      vals[5] = "total_s"

                    for (k = 1; k <= 5; k++) {
                        printf "%-12s", rows[k]
                        for (i = 1; i <= n; i++) {
                            ws = ws_sorted[i]
                            v = (vals[k] == "extract_s") ? extract_s[ws]
                              : (vals[k] == "sort_s")    ? sort_s[ws]
                              : (vals[k] == "merge_s")   ? merge_s[ws]
                              : (vals[k] == "align_s")   ? align_s[ws]
                              : total_s[ws]
                            printf "  %8.2fs ", v
                        }
                        printf "\n"
                    }
                }
            ' "$ROOT/pga-mpi-merge.timing.csv"
            echo
            echo "Raw CSV: $ROOT/pga-mpi-merge.timing.csv"
        fi
        ;;

    2)
        for cfg in "1:2" "4:8"; do
            per_node=${cfg%:*}; total=${cfg#*:}
            echo "=== 2 nodes × $per_node ranks/node (total=$total) ==="
            run_one "2x${per_node}" -N 2 --ntasks-per-node "$per_node" -n "$total"
        done
        ;;

    4)
        for cfg in "2:8" "4:16"; do
            per_node=${cfg%:*}; total=${cfg#*:}
            echo "=== 4 nodes × $per_node ranks/node (total=$total) ==="
            run_one "4x${per_node}" -N 4 --ntasks-per-node "$per_node" -n "$total"
        done
        ;;

    *)
        echo "ERROR: unsupported node count $NODES (script handles 1, 2, or 4)" >&2
        exit 1
        ;;
esac
