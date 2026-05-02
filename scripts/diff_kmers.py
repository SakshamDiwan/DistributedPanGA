#!/usr/bin/env python3
import sys
from collections import Counter

def load(path):
    with open(path) as f:
        return [line.rstrip("\n") for line in f]

def main():
    if len(sys.argv) != 3:
        print("Usage: diff_kmers.py <a.tuples> <b.tuples>", file=sys.stderr)
        sys.exit(2)
    a = load(sys.argv[1])
    b = load(sys.argv[2])
    sa, sb = set(a), set(b)
    only_a = sa - sb
    only_b = sb - sa
    print(f"A: {len(a)} lines ({len(sa)} unique)")
    print(f"B: {len(b)} lines ({len(sb)} unique)")
    print(f"Only in A: {len(only_a)}")
    print(f"Only in B: {len(only_b)}")
    if only_a or only_b:
        print("\nFirst 5 only in A:")
        for x in sorted(only_a)[:5]: print(f"  {x}")
        print("\nFirst 5 only in B:")
        for x in sorted(only_b)[:5]: print(f"  {x}")
        sys.exit(1)
    if a != b:
        ca, cb = Counter(a), Counter(b)
        diff_keys = [k for k in ca if ca[k] != cb[k]]
        print(f"\nMultiplicity differences: {len(diff_keys)}")
        for k in diff_keys[:5]:
            print(f"  {k}: A={ca[k]} B={cb[k]}")
        sys.exit(1)
    print("\nMATCH")

if __name__ == "__main__":
    main()
