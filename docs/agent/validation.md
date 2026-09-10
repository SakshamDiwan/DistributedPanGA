# Validation

Source commit inspected: `ce64406`. **No `make verify-*` target was executed in
this checkout** — the required fixtures do not exist here (see Prerequisites).
Every command below is source-derived from the Makefile and
`scripts/verify_*.sh` unless marked otherwise.

## Prerequisites

All five verifiers need `tests/fixtures/yeast/` artifacts, and the three MPI
ones additionally need an active Slurm allocation (they call `srun` and exit
non-zero if it is absent).

| Requirement | Needed by | Present here **[executed check]** |
|---|---|---|
| `tests/fixtures/yeast/yeast7.1gdb` | all five | **missing** (`tests/` does not exist) |
| `tests/fixtures/yeast/gold.tuples` | phase1, phase2, phase3-stage1, phase3-stage2 | **missing** |
| `tests/fixtures/yeast/yeast7.gix` | phase3-stage2-lcp | **missing** |
| `tests/fixtures/yeast/gold.lcp.kmers.simple` | phase3-stage2-lcp (else built inline) | **missing** |
| `tests/fixtures/yeast/gold.paf` | runtime script | **missing** |
| `~/FASTGA` tools (`FAtoGDB`, `GIXmake`, `GIXshow`, `FastGA`, `ALNtoPAF`) | gold generation, end-to-end | **missing** |
| active allocation / `srun` | phase3-stage1, phase3-stage2, phase3-stage2-lcp | not held in this session |

`tests/` is `.gitignore`d, so Git will never supply these. Regenerating them
means running `scripts/make_gold.sh`, which is a **deliberate, separate**
action — it deletes and rebuilds `yeast7.1gdb`/`yeast7.gix` in place and needs
`yeast7.fa.gz` plus `~/FASTGA`. Do not run it casually; preserve existing gold
artifacts.

## The five verifiers

Run from the repository root. Each `verify-*` target has the binary as a
prerequisite and the script rebuilds it too.

### `make verify-phase1`

Serial. Builds `pga-extract`, extracts from `yeast7`, external-`sort`s the
output, and requires it byte-identical to `gold.tuples`. No allocation needed.

### `make verify-phase2`

Serial. Builds `pga-sort`, which runs extract + `msd_sort`, then
external-`sort`s (because `msd_sort` orders by k-mer while `gold.tuples` is
lexicographic text order) and compares to `gold.tuples`. No allocation needed.
Note `pga-sort` calls `sort_records` only — its LCP byte is raw `msd_sort`
output, so this target says nothing about semantic LCPs.

### `make verify-phase3-stage1`

Needs an allocation. Builds `pga-mpi-stage1` and, for each `N` in
`RANK_COUNTS`, runs `srun -n $N` and compares rank 0's gathered output to
`gold.tuples` after an external sort.

### `make verify-phase3-stage2`

Needs an allocation. Builds `pga-mpi` and, for each `N`, runs `srun -n $N`,
concatenates `out.N$N.<r>.tuples` **in rank order**, then applies both checks —
both are mandatory:

- **Tier 2, order:** an `awk` walk over the *unsorted* rank-ordered
  concatenation asserting each line's k-mer (field `$4` in the default
  4-column dump) is `>=` the previous. This is what proves the distributed sort
  actually produced a globally ordered result rather than merely the right
  contents.
- **Tier 1, multiset:** `sort` the concatenation and `cmp` against
  `gold.tuples`. This proves nothing was lost or duplicated in the
  `Alltoallv`.

Neither check subsumes the other: Tier 1 passes for any permutation, Tier 2
passes for a sorted-but-lossy result. The script runs Tier 2 first and skips
Tier 1 for that `N` if it fails.

### `make verify-phase3-stage2-lcp`

Needs an allocation. **Keep this run separate from other builds**: the target
has no prerequisite and the script begins with
`make clean && make LCPS=1 pga-mpi`, deleting all six binaries. Do not
interleave it with a build or a run you care about.

It runs `PGA_DUMP_LCP=1 srun -n $N pga-mpi`, so the dump is the 5-column form
(`<contig> <position> <strand> <lcp> <hex_kmer>`; LCP is `$4`, k-mer is `$5`),
concatenates in rank order, and applies:

- **Check A** — for adjacent records with identical k-mers, byte 0 must be
  exactly `40`.
- **Check B** — for adjacent records with different k-mers, byte 0 must equal
  an independently implemented `simple_lcp` in `awk` (first differing base
  position, `[0,39]`).
- **Cross-check** — the derived `(lcp, distinct_kmer)` sequence must be
  byte-identical to `gold.lcp.kmers.simple`, or to an inline rebuild via
  `scripts/build_lcp_oracle_simple.sh` if the fixture is absent.

The oracle deliberately **does not** use GIXshow's `lcp` column, which is a
verbatim view of `msd_sort`'s internal byte 0 and disagrees with
bases-in-common at `part[]` boundaries (the script's header cites an observed
case: GIXshow shows 12 for `tttt...t -> caaa...a`, where bases-in-common is 0).
The oracle is tie-break independent because a transition LCP depends only on
the two distinct k-mers, not on ordering within either tie group.

On failure this script preserves its `mktemp` work directory and prints the
path; on success it removes it.

## `RANK_COUNTS`

The three MPI verifiers honour `RANK_COUNTS` (default `1 2 4 8`):

```bash
RANK_COUNTS="1 2" make verify-phase3-stage2
RANK_COUNTS=4     make verify-phase3-stage2-lcp
```

`RANK_COUNTS` has **no effect** on `scripts/runtime_test_phase4b.sh`.

Include `N=1` and at least one `N` that is not a power of two once fixtures
exist: the splitter's rounding path and the empty-slice paths are only
exercised at awkward rank counts.

## Recommended order inside one allocation

```bash
salloc --account=m4341 --constraint=cpu --qos=interactive --nodes=1 --time=00:30:00
cd /global/u2/h/hsc53/DistributedPanGA
make verify-phase1
make verify-phase2
make verify-phase3-stage1
make verify-phase3-stage2
make verify-phase3-stage2-lcp     # last: it runs `make clean` first
```

## `PGA_DUMP_LCP`

Read in `src/pga-mpi.c` only, and only the exact value `1` counts (the check is
`getenv(...)[0] == '1'`). Unset or any other value gives the 4-column dump.
Because the k-mer moves from `$4` to `$5`, a script written for one mode
silently misreads the other — `verify_phase3_stage2.sh` keys on `$4` and must
run *without* the variable; `verify_phase3_stage2_lcp.sh` sets it and keys on
`$5`.

## End-to-end runtime script (separate from the verifiers)

`scripts/runtime_test_phase4b.sh` is the only end-to-end driver. It is not a
`make` target and is not a correctness verifier in the sense above.

```bash
salloc --account=m4341 --constraint=cpu --qos=interactive --nodes=1 --time=00:30:00
bash scripts/runtime_test_phase4b.sh      # NOT under srun
```

It builds `pga-mpi-merge`, sweeps world sizes chosen from
`$SLURM_JOB_NUM_NODES` (1 → `-n 1,2,4,8`; 2 → `2x1` and `2x4`; 4 → `4x2` and
`4x4`; anything else is an error), compares `sort`ed concatenated
`out.rank*.paf` against `sort`ed `gold.paf` per configuration, and prints a
stage × world-size table from `pga-mpi-merge.timing.csv`.

Caveats to keep in mind when reading its output — see `perlmutter.md` for the
full list: it truncates the CSV and deletes `out.rank*.paf` on every
invocation, uses fixed `/tmp/gold.sorted.paf` and `/tmp/got.sorted.paf`, runs
from the repository root, and requires `/global/u2/s/sdiwan/FASTGA` to exist.
Its `TOTAL` row is the sum of per-stage `MPI_MAX` values, not measured wall
time.

## Gold artifact provenance

Record all of this alongside any result derived from gold data.

| Artifact | Produced by | Notes |
|---|---|---|
| `yeast7.1gdb` | `FAtoGDB yeast7.fa.gz` | input is `yeast7.fa.gz`, not in Git |
| `yeast7.gix` | `GIXmake yeast7` | rebuilt again by step 4 of `make_gold.sh` |
| `gold.kmers` | `GIXshow yeast7.gix \| sort` | |
| `gold.lcp.kmers.simple` | `scripts/build_lcp_oracle_simple.sh yeast7.gix` | GIXshow k-mers, LCPs recomputed as bases-in-common; uses `$GIXSHOW`, default `~/FASTGA/GIXshow` |
| `gold.paf` | `FastGA -1:gold yeast7` then `ALNtoPAF gold.1aln` | **one** genome argument = self-alignment mode (`SELF=1`); `yeast7 yeast7` would be cross-mode and give roughly half the alignments |
| `gold.tuples` | `notes/GIXmake_dump -v yeast7` with `GIXMAKE_DUMP_PATH=/tmp/yeast_dump.tuples`, then `sort` | instrumented GIXmake; step 5 asserts byte-identity against `GIXshow \| normalize_gixshow \| sort` |

`notes/GIXmake_dump` is a **tracked prebuilt binary**; rebuilding it via
`notes/build_dump.sh` requires `~/FASTGA` sources. Treat `notes/` as reference
and instrumentation, not production source.

Parameters baked into the gold set: `KMER=40`, `TMER=12`, `SMER=8`, and
FastGA's own defaults for the alignment stage. Any change to those constants
invalidates every gold artifact.

## Reporting results honestly

- Say which target ran, at which `RANK_COUNTS`, and paste the PASS/FAIL lines.
- Do not describe a source-derived command as executed.
- The historical yeast7 figure of **43,126 alignments** belongs to a specific
  input, parameter set and reference. It has **not** been reproduced in this
  checkout and must not be cited as if it had.
- A zero exit status from `pga-mpi-merge` does not imply PAF output — check the
  files.
