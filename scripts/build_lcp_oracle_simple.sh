#!/bin/bash
# Build LCP oracle using simple base-level LCP calculation.
# Extracts distinct kmers from GIXshow output and recalculates LCPs
# using the same simple_kmer_lcp_bases logic as our MPI pipeline.
#
# Usage:  build_lcp_oracle_simple.sh <gix-file> <output-file>

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
KBYTES=10

# Extract k-mers from GIXshow and recalculate LCPs with simple base-level logic
"$GIXSHOW" "$GIX_FILE" | gawk -v KMER=$KMER -v KBYTES=$KBYTES '
    BEGIN {
        # base->2bit map
        b2v["a"] = 0; b2v["c"] = 1; b2v["g"] = 2; b2v["t"] = 3
        hex = "0123456789abcdef"
    }
    
    # Convert DNA string to byte array (10 bytes)
    function dna_to_bytes(s, bytes,    i, k, byte) {
        for (i = 0; i < KBYTES; i++) {
            byte = 0
            for (k = 0; k < 4; k++)
                byte = byte * 4 + b2v[substr(s, i*4 + k + 1, 1)]
            bytes[i] = byte
        }
    }
    
    # Convert byte array to hex string
    function bytes_to_hex(bytes,    out, i) {
        out = ""
        for (i = 0; i < KBYTES; i++)
            out = out substr(hex, int(bytes[i]/16) + 1, 1) substr(hex, (bytes[i] % 16) + 1, 1)
        return out
    }
    
    # Simple base-level LCP: count matching bases (0-40)
    function simple_lcp(bytes1, bytes2,    byte_idx, base_in_byte, diff, shift) {
        for (byte_idx = 0; byte_idx < KBYTES; byte_idx++) {
            if (bytes1[byte_idx] != bytes2[byte_idx]) {
                diff = xor(bytes1[byte_idx], bytes2[byte_idx])
                for (base_in_byte = 0; base_in_byte < 4; base_in_byte++) {
                    shift = (3 - base_in_byte) * 2
                    if (and(rshift(diff, shift), 3)) {
                        return byte_idx * 4 + base_in_byte
                    }
                }
            }
        }
        return KMER  # All 40 bases match
    }
    
    NR == 1 { next }   # skip header
    
    {
        # Parse GIXshow line: "  Index: kmer mask lcp sign contig | position"
        idx_pos = index($0, ":")
        if (idx_pos == 0) next
        rest = substr($0, idx_pos + 1)
        n = split(rest, f, /[ \t]+/)
        kmer = (f[1] == "") ? f[2] : f[1]
        if (length(kmer) != KMER) next
        
        # Convert kmer to bytes for comparison
        dna_to_bytes(kmer, curr_bytes)
        
        if (NR == 2) {
            # First distinct kmer - LCP is always 0
            prev_hex = bytes_to_hex(curr_bytes)
            print 0, prev_hex
            for (i = 0; i < KBYTES; i++) prev_bytes[i] = curr_bytes[i]
            next
        }
        
        # Check if same as previous (GIXshow uses "*" for identical)
        is_same = 1
        for (i = 0; i < KBYTES; i++) {
            if (curr_bytes[i] != prev_bytes[i]) {
                is_same = 0
                break
            }
        }
        
        if (!is_same) {
            # Different kmer - calculate LCP using simple logic
            lcp = simple_lcp(prev_bytes, curr_bytes)
            hex_kmer = bytes_to_hex(curr_bytes)
            print lcp, hex_kmer
            for (i = 0; i < KBYTES; i++) prev_bytes[i] = curr_bytes[i]
        }
        # Skip identical kmers (would have LCP=40 in our encoding)
    }
' > "$OUTPUT"

echo "Wrote $(wc -l < "$OUTPUT") distinct-kmer LCPs (simple encoding) to $OUTPUT"