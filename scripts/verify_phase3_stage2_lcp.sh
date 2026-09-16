#!/bin/bash
# Phase 3 stage 2 LCP regression: build pga-mpi with LCPS=1, run at
# N ∈ {1,2,4,8}, concatenate per-rank outputs in rank order, and run two
# independent LCP checks plus a cross-check against a simple-LCP oracle.
#
# Encoding used throughout this test: byte 0 = number of matching bases with
# the previous record's k-mer, range [0, 40]. Identical k-mers -> 40,
# distinct k-mers -> 0..39 based on the position of the first differing base.
# This is a *semantic* LCP (true bases-in-common), NOT msd_sort's internal
# byte-0 encoding. Both pga-mpi (via recalc_all_lcps) and the oracle
# (build_lcp_oracle_simple.sh) use this same convention.
#
# We deliberately do NOT compare against GIXshow's "lcp" column, because that
# column is a verbatim display of msd_sort's internal byte 0 — which has
# implementation-specific encodings at part[] boundaries that don't reflect
# semantic bases-in-common (verified empirically: GIXshow shows 12 for the
# transition tttt...t -> caaa...a where actual bases-in-common is 0). See
# implementation_plan.md §7a.1b for the full encoding analysis.
#
# Check A — within-group invariant:
#   For every adjacent pair where kmer[i] == kmer[i-1], byte 0 must equal 40
#   (all 40 bases match).
#
# Check B — between-group LCP:
#   For every adjacent pair where kmer[i] != kmer[i-1], byte 0 must equal
#   simple_lcp(prev, curr) — the position of the first differing base in [0, 39].
#
# Cross-check — simple-LCP oracle:
#   gold.lcp.kmers.simple is built by build_lcp_oracle_simple.sh, which extracts
#   distinct k-mers from GIXshow and recomputes LCPs with the same simple_lcp
#   logic. Our concatenated per-rank output, reduced to (distinct_kmer, lcp)
#   pairs in sort order, must cmp byte-identical to that oracle.
#
# REQUIRES: an active Slurm allocation. Uses srun.
# Exits 0 on success, non-zero on any failure.

set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
# Overridable so the suite can run against a small synthetic fixture as well
# as the yeast reference. Defaults are the historical values.
FIXTURES=${PGA_FIXTURES:-$ROOT/tests/fixtures/yeast}
GDB_STEM=${PGA_GDB_STEM:-$FIXTURES/yeast7}
GIX_FILE=${PGA_GIX_FILE:-$GDB_STEM.gix}
ORACLE_FIXTURE=$FIXTURES/gold.lcp.kmers.simple
WORK=$(mktemp -d)
# Preserve $WORK on failure for postmortem; clean up on full success only.
trap '[ "$overall_status" -eq 0 ] && rm -rf "$WORK" || echo "[work dir preserved at $WORK for postmortem]" >&2' EXIT
overall_status=0

KMER=40

if [ ! -f "$GDB_STEM.1gdb" ]; then
    echo "ERROR: $GDB_STEM.1gdb not found." >&2
    exit 1
fi
if [ ! -f "$GIX_FILE" ]; then
    echo "ERROR: $GIX_FILE not found. Run scripts/make_gold.sh first." >&2
    exit 1
fi
if ! command -v srun >/dev/null 2>&1; then
    echo "ERROR: srun not found. Run inside a Slurm allocation." >&2
    exit 1
fi

echo "[build] pga-mpi (LCPS=1 is the Makefile default; passed explicitly for clarity)..."
( cd "$ROOT" && make clean >/dev/null && make LCPS=1 pga-mpi >/dev/null )

# -----------------------------------------------------------------------------
# Locate the simple-LCP oracle. Preferred path:
#   tests/fixtures/yeast/gold.lcp.kmers.simple
# (persisted artifact from scripts/make_gold.sh). Fallback: build inline into
# a tempfile if the fixture is missing — keeps the script self-contained for
# fresh checkouts that haven't regenerated gold yet.
#
# The oracle uses simple_lcp encoding (bases-in-common, [0, 40]), NOT GIXshow's
# raw byte 0. See header comment for why.
# -----------------------------------------------------------------------------
if [ -f "$ORACLE_FIXTURE" ]; then
    echo "[oracle] using cached fixture $ORACLE_FIXTURE"
    ORACLE=$ORACLE_FIXTURE
else
    echo "[oracle] $ORACLE_FIXTURE missing; building inline (consider running scripts/make_gold.sh)"
    ORACLE="$WORK/oracle.lcp.txt"
    "$ROOT/scripts/build_lcp_oracle_simple.sh" "$GIX_FILE" "$ORACLE"
fi
echo "  oracle has $(wc -l < "$ORACLE") distinct k-mer transitions"

# -----------------------------------------------------------------------------
# AWK program shared by all rank counts. Reads the 5-column LCP-aware dump,
# performs Check A and Check B, and emits a (distinct_kmer, first_lcp) sequence
# to stdout for cross-checking against the oracle.
#
# Input format:  <contig> <position> <strand> <lcp> <hex_kmer>
# (see dump_records_as_tuples in src/pga-mpi.c with PGA_DUMP_LCP=1)
#
# Uses simple base-level LCP calculation that counts actual matching bases (0-40),
# matching our recalc_all_lcps implementation.
# -----------------------------------------------------------------------------
read -r -d '' CHECKER_AWK << 'AWK_EOF' || true
BEGIN {
    KMER = 40
    KBYTES = 10

    # Hex digit lookup.
    hexstr = "0123456789abcdef"
    for (i = 0; i < 16; i++) hexval[substr(hexstr, i+1, 1)] = i

    a_violations = 0
    b_violations = 0
    a_first_bad = -1
    b_first_bad = -1
    line_num = 0
    distinct_count = 0
}

# Parse a 20-char hex string into 10 byte values; store in array kbyte[0..9].
function parse_kmer(hex,    i, hi, lo) {
    for (i = 0; i < KBYTES; i++) {
        hi = hexval[substr(hex, 2*i+1, 1)]
        lo = hexval[substr(hex, 2*i+2, 1)]
        kbyte[i] = hi * 16 + lo
    }
}

# Calculate simple base-level LCP: number of matching bases (0-40)
# between k-mers stored in arrays a[] and b[].
function simple_lcp(    byte_idx, base_in_byte, diff, shift) {
    for (byte_idx = 0; byte_idx < KBYTES; byte_idx++) {
        if (a[byte_idx] != b[byte_idx]) {
            diff = xor(a[byte_idx], b[byte_idx])
            for (base_in_byte = 0; base_in_byte < 4; base_in_byte++) {
                shift = (3 - base_in_byte) * 2
                if (and(rshift(diff, shift), 3)) {
                    return byte_idx * 4 + base_in_byte
                }
            }
        }
    }
    return 40  # All 40 bases match
}

{
    line_num++
    lcp = $4 + 0
    hex = $5
    parse_kmer(hex)

    if (line_num == 1) {
        # First record: LCP should be 0 (no predecessor).
        if (lcp != 0) {
            a_violations++
            if (a_first_bad < 0) a_first_bad = line_num
            printf("Check A FAIL line 1: first record LCP=%d, expected 0 (kmer=%s)\n",
                   lcp, hex) > "/dev/stderr"
        }
        for (i = 0; i < KBYTES; i++) prev[i] = kbyte[i]
        # Emit oracle line for first record
        distinct_count++
        printf("%d %s\n", 0, hex)
        prev_hex = hex
        next
    }

    # Compare current to previous.
    same = 1
    for (i = 0; i < KBYTES; i++) {
        if (kbyte[i] != prev[i]) { same = 0; break }
    }

    if (same) {
        # Check A: within-group byte 0 must be 40 (identical k-mers).
        if (lcp != 40) {
            a_violations++
            if (a_first_bad < 0) a_first_bad = line_num
            if (a_violations <= 5)
                printf("Check A FAIL line %d: identical kmers LCP=%d, expected 40 (kmer=%s)\n",
                       line_num, lcp, hex) > "/dev/stderr"
        }
    } else {
        # Check B: between-group LCP must match our simple calculator.
        for (i = 0; i < KBYTES; i++) a[i] = prev[i]
        for (i = 0; i < KBYTES; i++) b[i] = kbyte[i]
        expected = simple_lcp()
        if (lcp != expected) {
            b_violations++
            if (b_first_bad < 0) b_first_bad = line_num
            if (b_violations <= 5)
                printf("Check B FAIL line %d: LCP=%d, expected %d (prev=%s curr=%s)\n",
                       line_num, lcp, expected, prev_hex, hex) > "/dev/stderr"
        }
        # Emit for oracle cross-check (simple-LCP encoding, [0, 40] range,
        # same as build_lcp_oracle_simple.sh produces).
        distinct_count++
        printf("%d %s\n", lcp, hex)
        for (i = 0; i < KBYTES; i++) prev[i] = kbyte[i]
        prev_hex = hex
    }
}

END {
    printf("Check A: %d within-group violations (first at line %d)\n",
           a_violations, a_first_bad) > "/dev/stderr"
    printf("Check B: %d between-group violations (first at line %d)\n",
           b_violations, b_first_bad) > "/dev/stderr"
    printf("Distinct k-mer count: %d\n", distinct_count) > "/dev/stderr"
    if (a_violations > 0 || b_violations > 0) exit 1
}
AWK_EOF

# -----------------------------------------------------------------------------
# Per-rank-count loop.
# -----------------------------------------------------------------------------
RANK_COUNTS="${RANK_COUNTS:-1 2 4 8}"

for N in $RANK_COUNTS; do
    echo "=== N=$N ==="
    OUT_STEM="$WORK/out.N$N"
    COMBINED="$WORK/combined.N$N.tuples"
    OUR_ORACLE="$WORK/our_oracle.N$N.txt"

    echo "  [run] PGA_DUMP_LCP=1 srun -n $N pga-mpi ..."
    PGA_DUMP_LCP=1 srun -n "$N" "$ROOT/pga-mpi" "$GDB_STEM" "$OUT_STEM"

    echo "  [cat] concatenating ${N} per-rank outputs in rank order..."
    : > "$COMBINED"
    for ((r = 0; r < N; r++)); do
        cat "$OUT_STEM.$r.tuples" >> "$COMBINED"
    done

    echo "  [check A+B] running awk checker..."
    if gawk "$CHECKER_AWK" "$COMBINED" > "$OUR_ORACLE" 2> "$WORK/check.err"; then
        echo "  Check A+B PASS"
        cat "$WORK/check.err"
    else
        echo "  Check A or B FAIL" >&2
        cat "$WORK/check.err" >&2
        overall_status=1
        continue
    fi

    echo "  [cross-check] vs simple-LCP oracle..."
    if cmp -s "$OUR_ORACLE" "$ORACLE"; then
        echo "  Cross-check PASS: $(wc -l < "$OUR_ORACLE") distinct-kmer LCPs match oracle"
    else
        echo "  Cross-check FAIL: our LCPs differ from oracle" >&2
        echo "    ours:   $(wc -l < "$OUR_ORACLE") lines" >&2
        echo "    oracle: $(wc -l < "$ORACLE") lines" >&2
        diff <(head -10 "$OUR_ORACLE") <(head -10 "$ORACLE") | head -20 >&2 || true
        overall_status=1
    fi
done

if [ "$overall_status" -eq 0 ]; then
    echo ""
    echo "ALL PASS: stage 2 LCPs verified at N ∈ { $RANK_COUNTS } (Check A + Check B + simple-LCP cross-check)"
else
    echo ""
    echo "SOME FAILED" >&2
fi
exit "$overall_status"