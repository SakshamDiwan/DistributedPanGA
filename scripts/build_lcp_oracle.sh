#!/bin/bash
# Build the LCP oracle: a sequence of "<lcp> <hex_kmer>" lines, one per
# distinct k-mer in GIXshow output, in GIXshow order. Used by the Phase 3
# stage 2 LCP verifier. See scripts/verify_phase3_stage2_lcp.sh for the
# checks this oracle backs.
#
# Usage:  build_lcp_oracle.sh <gix-file> <output-file>
#
# Why this is tie-break-independent:
#   GIXshow's LCP column gives the byte-resolution LCP between adjacent
#   records in sort order. Within a tie group of identical k-mers, every
#   such LCP is KMER (40), printed as "*". At a transition from one
#   distinct k-mer to a different one, the LCP value depends only on the
#   two k-mers, NOT on the tie order within either group. So the sequence
#   of (transition_lcp, new_kmer) pairs is fully determined by the set of
#   distinct k-mers and their sort order — invariant under any tie-break
#   permutation that any sorter (FastGA's or ours) might choose.

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <gix-file> <output-file>" >&2
    exit 1
fi

GIX_FILE=$1
OUTPUT=$2
GIXSHOW=${GIXSHOW:-$HOME/FASTGA/GIXshow}

if [ ! -f "$GIX_FILE" ]; then
    echo "ERROR: $GIX_FILE not found." >&2
    exit 1
fi
if [ ! -x "$GIXSHOW" ]; then
    echo "ERROR: $GIXSHOW not found or not executable." >&2
    exit 1
fi

KMER=40

"$GIXSHOW" "$GIX_FILE" | gawk -v KMER=$KMER '
    BEGIN {
        # base->2bit map (a/c/g/t -> 0/1/2/3, lowercase only since GIXshow always
        # emits lowercase).
        b2v["a"] = 0; b2v["c"] = 1; b2v["g"] = 2; b2v["t"] = 3
        hex = "0123456789abcdef"
    }
    # Pack a 40-base DNA string into 20 hex chars (10 bytes, 4 bases each).
    function pack_kmer(s,    out, i, k, byte) {
        out = ""
        for (i = 0; i < KBYTES; i++) {
            byte = 0
            for (k = 0; k < 4; k++)
                byte = byte * 4 + b2v[substr(s, i*4 + k + 1, 1)]
            out = out substr(hex, int(byte/16) + 1, 1) substr(hex, (byte % 16) + 1, 1)
        }
        return out
    }
    BEGIN { KBYTES = KMER / 4 }

    NR == 1 { next }   # skip header
    {
        # GIXshow line: "  Index: kmer mask lcp sign contig | position"
        idx_pos = index($0, ":")
        if (idx_pos == 0) next
        rest = substr($0, idx_pos + 1)
        n = split(rest, f, /[ \t]+/)
        # f[1] is empty (leading space split), f[2] is kmer, f[3] is mask, f[4] is lcp.
        kmer = (f[1] == "") ? f[2] : f[1]
        lcp_field = (f[1] == "") ? f[4] : f[3]
        if (length(kmer) != KMER) next

        if (NR == 2) {
            prev_kmer = kmer
            print 0, pack_kmer(kmer)
            next
        }
        if (kmer != prev_kmer) {
            # Transition: lcp_field is a number 1..KMER-1.
            print lcp_field+0, pack_kmer(kmer)
            prev_kmer = kmer
        }
        # Within-group lines (lcp_field = "*") contribute nothing.
    }
' > "$OUTPUT"

echo "Wrote $(wc -l < "$OUTPUT") distinct-kmer LCPs to $OUTPUT"
