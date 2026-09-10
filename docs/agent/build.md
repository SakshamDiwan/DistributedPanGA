# Build

Source commit inspected: `ce64406`. Verification status is stated per command:
**[executed]** means it was run in this checkout on a Perlmutter login node on
2026-09-10; **[source-derived]** means it was read out of the Makefile or a
script and not run.

## Targets

`make all` builds six binaries **[executed, all six produced]**:

| Target | Compiler | Sources beyond core |
|---|---|---|
| `pga-extract` | literal `gcc` | — |
| `pga-sort` | literal `gcc` | `sort.c` |
| `pga-mpi-stage1` | `$(CC)` | `sort.c dsort.c` |
| `pga-mpi` | `$(CC)` | `sort.c dsort.c` |
| `pga-merge` | literal `gcc` | `sort.c kmer_adapter.c fastga_pipeline.c` |
| `pga-mpi-merge` | `$(CC)` | `sort.c dsort.c kmer_adapter.c fastga_pipeline.c contig_assignment.c` |

Core (`CORE_SRCS`) is `extract.c pack.c tables.c work_split.c record.c`.
Every target also compiles all of `LIB_SRCS`: `GDB.c gene_core.c libfastk.c
ONElib.c MSDsort.c align.c alncode.c ANO.c RSDsort.c`.

`pga-mpi-merge` is the final executable. `pga-mpi` is the distributed-sorting
diagnostic. There are no object files and no static library — every target
recompiles all sources into one binary, so builds are whole-program and take
roughly a minute each.

`make clean` removes all six binaries. It does **not** touch outputs
(`out*.paf`, `pga-mpi-merge.timing.csv`, `logs/`).

## Flags

```
CC     = mpicc
CFLAGS = -O3 -Wall -Wno-unused-result -fno-strict-aliasing -I./lib -I./src
LIBS   = -lpthread -lm -lz
LCPS  ?= 1        # LCPS=1 appends -DLCPs
```

- `-fno-strict-aliasing` is required: the adapter structs are cast to FastGA's
  `Kmer_Stream`/`Post_List` and seed payloads are punned through `uint8*`.
- `-DLCPs` is the default and must stay on for anything compared against
  FastGA. It selects semantic-LCP recalculation in `recalc_all_lcps` and the
  cross-rank seam fixup; without it, byte 0 becomes a 0/1 boundary marker.
- Because `CC` is only referenced by the three MPI recipes, **`CC=...` changes
  nothing for `pga-extract`, `pga-sort`, or `pga-merge`** — those hardcode
  `gcc`.

## Environment on Perlmutter

Loaded modules in this checkout's shell **[executed: `module list`]**:
`PrgEnv-gnu/8.7.0`, `gcc-native/14`, `cray-mpich/9.1.0`, `cray-libsci/26.03.0`,
`craype/2.7.36`, `craype-x86-milan`, plus `cudatoolkit/13.2`,
`craype-accel-nvidia80` and `gpu/1.0` (the default login-node GPU set — this is
a CPU-only project, so those are inert here).

Versions observed **[executed]**:

| Tool | Result |
|---|---|
| `gcc --version` | `gcc (SUSE Linux) 14.3.0` |
| `cc --version` | `gcc-14 (SUSE Linux) 14.3.0` (Cray wrapper) |
| `mpicc -show` | wraps `gcc` against `/opt/cray/pe/mpich/**9.0.1**/ofi/gnu/**12.3**` |

Note the mismatch: `mpicc` on `PATH` comes from the MPICH **9.0.1 / gnu-12.3**
tree while the loaded module is `cray-mpich/9.1.0`. Both build and link
successfully, but only the Cray wrapper `cc` is consistent with the loaded
modules, so prefer:

```bash
make CC=cc pga-mpi-merge
```

**[executed: `make -B CC=cc pga-mpi` completed with 0 errors and produced a
409128-byte binary versus 388888 bytes from `mpicc`]** — the size difference
confirms a different MPI library is being linked. Use one or the other
consistently; do not mix binaries from the two wrappers in one measurement.

Dependencies are all satisfied by the default environment: pthreads, libm,
zlib. Nothing needs to be installed to build.

## Rebuild pitfalls

Recipes depend only on `.c` files. **Headers are not prerequisites, and neither
are `LCPS` or `CFLAGS`.** Editing `src/constants.h`, `src/record.h`, or
changing a flag therefore leaves a stale binary that `make` calls up to date.

Verified live in this checkout **[executed]**: with a fresh `pga-mpi` present,
`make -n LCPS=0 pga-mpi` and `make -n CFLAGS=-O2 pga-mpi` both printed
`make: 'pga-mpi' is up to date.` — no recompile, silently wrong binary.

Always force the rebuild when a header or flag changed:

```bash
make -B LCPS=1 pga-mpi-merge
```

Second pitfall: a command-line `CFLAGS` **replaces** the Makefile's value, and
because `-DLCPs` is appended to `CFLAGS` it is lost with it. Verified
**[executed]**:

```
$ make -Bn CFLAGS=-O2 pga-mpi | head -1
mpicc -O2 -o pga-mpi src/pga-mpi.c ...            # no -DLCPs, no -I, no -fno-strict-aliasing
```

So `make CFLAGS=-O2 LCPS=1 ...` still produces an **LCPS=0-equivalent** binary
that will fail the LCP checks, and drops `-I./lib -I./src` and
`-fno-strict-aliasing` as well. Always inspect the expanded command with
`make -Bn <target>` before trusting a flag override. To add flags rather than
replace them, edit the Makefile or pass them through a target-specific
variable — not a bare `CFLAGS=`.

For a genuinely clean state:

```bash
make clean && make -B CC=cc all
```

## Expected warnings

The build is warning-noisy and that is pre-existing, not a regression. From
`lib/` (vendored FastGA): `-Wformat-overflow`, `-Wrestrict`,
`-Wmaybe-uninitialized`, `-Wstringop-overflow` in `align.c` and `ANO.c`.

From our own sources, compiling `src/pga-mpi-merge.c` alone **[executed]**
yields exactly two warnings, both dead code:

```
src/pga-mpi-merge.c:555:1: warning: label 'cleanup' defined but not used
src/pga-mpi-merge.c:102:14: warning: 'read_file_all' defined but not used
```

Neither affects behaviour; both are listed as follow-ups in `state.md`. Treat
any *new* warning from `src/` as a real finding.

## Outputs and where they land

Build products are the six binaries in the repository root; all are
`.gitignore`d. Run products (`out*.paf`, `pga-mpi-merge.timing.csv`, `logs/`)
are also ignored and land in the **current working directory**, not next to the
binary — see `perlmutter.md` before running from a shared filesystem.

## Verification summary

| Command | Status |
|---|---|
| `make all` | **executed**, exit 0, six binaries produced |
| `make -B CC=cc pga-mpi` | **executed**, exit 0 |
| `make -n {LCPS=0,CFLAGS=-O2} pga-mpi` | **executed**, both reported "up to date" |
| `make -Bn CFLAGS=-O2 pga-mpi` | **executed**, `-DLCPs` confirmed absent |
| `make verify-phase*` | **not executed** — fixtures absent, see `validation.md` |
| any `srun` run | **not executed** — no allocation held in this session |
