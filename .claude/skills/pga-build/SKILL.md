---
name: pga-build
description: Build DistributedPanGA targets on Perlmutter and diagnose stale-binary or flag-override problems. Use when asked to build, rebuild, or clean pga-extract, pga-sort, pga-mpi-stage1, pga-mpi, pga-merge, or pga-mpi-merge, or when a binary seems not to reflect a source or flag change.
---

Full reference, including verified compiler/MPI versions and both rebuild
pitfalls: `docs/agent/build.md`. Correctness constraints: `AGENTS.md`.

## Workflow

1. Identify the target. `pga-mpi-merge` is the final end-to-end executable;
   `pga-mpi` is the distributed-sorting diagnostic. `make all` builds all six.
2. Prefer the Cray wrapper, which matches the loaded `cray-mpich` module:
   `make CC=cc <target>`. Note `CC` affects only the three MPI targets —
   `pga-extract`, `pga-sort` and `pga-merge` hardcode `gcc`.
3. If a header (`src/*.h`), `LCPS`, or `CFLAGS` changed, force the rebuild:
   `make -B LCPS=1 <target>`. Headers and flags are **not** Makefile
   prerequisites, so plain `make` will report "up to date" and leave a stale
   binary.
4. Before trusting any flag override, print the expanded command with
   `make -Bn <target>` and confirm `-DLCPs`, `-I./lib -I./src` and
   `-fno-strict-aliasing` are all present. A command-line `CFLAGS=` **replaces**
   the Makefile value and silently drops all three.
5. Ignore warnings from `lib/` — they are pre-existing upstream FastGA noise.
   Treat any new warning from `src/` as a real finding.
6. For a clean, benchmark-ready state: `make clean && make -B CC=cc all`. Never
   mix `cc`- and `mpicc`-built binaries in one measurement; they link different
   MPI libraries.
