# PR 2 results — large-count exchange (R-06)

- Job 59713118 (2026-10-12, nid004429), commit `fd6e104` (before the rebase
  onto `da9830d`/`d470c6b`; the code is the same apart from PR 1's
  measurement follow-up).
- PR 1 reference binary built from `pr1/align-instrumentation` in the same
  job.
- Raw data: `$SCRATCH/pga-pr2-59713118/`.

## Correctness (all green)

| Check | Result |
|---|---|
| `make test` (44,231 checks incl. R-06 1,733), `make test-asan` incl. R-06 | pass (login node) |
| `test-mpi` ws 1/2/3/4/8, default budget | pass, 0 FAIL lines |
| `test-mpi` ws 1/2/3/4/8, `PGA_EXCHANGE_ROUND_BYTES=256` (many stage-2 rounds) | pass, 0 FAIL lines |
| small6 e2e N=1/2/4/8, default and 4096-byte budget | ALL PASS both |
| yeast7 `verify-phase3-stage2` (Tier 1 + Tier 2), 1 MiB budget, N=1/2/4/8 | ALL PASS |
| yeast7 N=1/2/8 and scerevisiae8 N=4: PR 2 at 8 MiB budget vs PR 1 default | gold MATCH; **raw per-rank PAF byte-identical** in all 4 |

The byte-identical raw PAFs are the "exact exchange semantics" evidence. They
pass even with many rounds in both exchanges, and they show the chunked path
writes exactly what one unbounded `Alltoallv` wrote.

## Rounds used at an 8 MiB budget (measured)

| run | stage-2 rounds | seed N / C rounds | seed staging |
|---|---|---|---|
| yeast7 N=1 | 123 | — (single-rank shortcut) | 0 |
| yeast7 N=2 | 63 | 169 / 10 | 16 MiB |
| yeast7 N=8 | 18 | 48 / 3 | 16 MiB |
| scerevisiae8 N=4 | 39 | 141 / 7 | 16 MiB |

Wall time is unchanged within noise despite the rounds (yeast7 N=1 64.8 vs
65.5 s, N=8 10.52 vs 10.56 s). At the default 1 GiB budget, yeast7 and
scerevisiae8 at N ≥ 4 take the single-`Alltoallv` path for stage 2.

## Memory (measured, `[pr1-mem]` per-stage peaks, worst rank)

| | PR 1 seed-redistribution peak | PR 2 seed-redistribution peak | PR 2 overall peak (stage) |
|---|---|---|---|
| yeast7 N=2 | 2,526 MB | 40 MB | 1,551 MB (stage 2) |
| yeast7 N=8 | 628 MB | 47 MB | 412 MB (stage 2) |
| scerevisiae8 N=4 | 2,236 MB | ~50 MB | 1,080 MB (alignment) |

Streaming the seed files removes the full send and receive buffers. The
overall peak per rank drops another ~35–50% at N>1 on top of PR 1. With the
redistribution buffers gone, stage 2 or the alignment sort array now sets the
peak.

## What this unblocks

scerevisiae8 at N=2 aborted on `main` and PR 1 (seed exchange > `INT_MAX`).
The primate datasets need this at every planned rank count (Step 3).
