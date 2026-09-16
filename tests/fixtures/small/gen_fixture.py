#!/usr/bin/env python3
"""Generate a small multi-contig FASTA fixture for the end-to-end tests.

The fixture must produce a MEANINGFUL number of alignments -- a random
sequence self-aligns to almost nothing, so every assertion downstream would
pass vacuously. So we plant two families of near-repeats across different
subsets of contigs, at different offsets and with a few percent divergence.
That gives FastGA real, non-trivial work: cross-contig alignments that are
similar but not identical.

Deterministic: the seed is fixed, so the fixture and its gold artifacts are
reproducible. Output is ~48 kb over 6 contigs, small enough to run in seconds.

Usage:  gen_fixture.py <out.fa> [seed]
"""
import random
import sys

BASES = "ACGT"
NCONTIG = 6
CONTIG_LEN = 8000
FAMILY_A_LEN = 1500       # planted in contigs 0, 2, 4
FAMILY_B_LEN = 1200       # planted in contigs 1, 3, 5
DIVERGENCE = 0.02         # per-base substitution rate between copies


def rand_seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def mutate(rng, seq, rate):
    out = []
    for c in seq:
        if rng.random() < rate:
            out.append(rng.choice([b for b in BASES if b != c]))
        else:
            out.append(c)
    return "".join(out)


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: gen_fixture.py <out.fa> [seed]")
    out_path = sys.argv[1]
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20260912
    rng = random.Random(seed)

    family_a = rand_seq(rng, FAMILY_A_LEN)
    family_b = rand_seq(rng, FAMILY_B_LEN)

    contigs = []
    for i in range(NCONTIG):
        seq = list(rand_seq(rng, CONTIG_LEN))
        # plant one family copy per contig, alternating, at a varying offset
        fam = family_a if i % 2 == 0 else family_b
        copy = mutate(rng, fam, DIVERGENCE)
        off = 500 + (i * 900) % (CONTIG_LEN - len(copy) - 600)
        seq[off:off + len(copy)] = list(copy)
        # a second, shorter tandem copy inside the contig, for within-contig hits
        short = mutate(rng, fam[:600], DIVERGENCE)
        off2 = off + len(copy) + 200
        if off2 + len(short) < CONTIG_LEN:
            seq[off2:off2 + len(short)] = list(short)
        contigs.append("".join(seq))

    with open(out_path, "w") as f:
        for i, s in enumerate(contigs):
            f.write(">contig_%d\n" % i)
            for j in range(0, len(s), 80):
                f.write(s[j:j + 80] + "\n")

    total = sum(len(s) for s in contigs)
    print("wrote %s: %d contigs, %d bp, seed %d" % (out_path, len(contigs), total, seed))


if __name__ == "__main__":
    main()
