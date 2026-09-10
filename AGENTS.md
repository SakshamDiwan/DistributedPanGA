# DistributedPanGA — agent instructions

MPI distributification of FastGA's pangenome self-alignment. Each rank extracts
syncmer-selected 40-mer seed records from a shared GDB, redistributes them by
k-mer prefix, sorts locally, then merges/chains/aligns its assigned contig
pairs. Final executable: `pga-mpi-merge`. Inspected source commit: `ce64406`.

## Environment

- NERSC Perlmutter. CPU account **m4341** (`m4341_g` is the GPU account; GPU
  work is out of scope — see `docs/agent/research-plan.md`).
- Interactive start: `salloc --account=m4341 --constraint=cpu --qos=interactive --nodes=1 --time=00:30:00`
- CPU node: 128 physical cores / 256 logical, 512 GB. `make all` builds on a
  login node in ~1 min.
- No `tests/` fixtures, no `~/FASTGA`, and no `ALNtoPAF`/`FastGA` on `PATH` in
  this checkout. Anything needing gold data or the FastGA tools cannot run yet.

## Correctness constraints (do not silently change)

- `KMER=40`, `KBYTES=10`, `TMER=12`, `SMER=8`, `NUM_BUCK=1024` (`src/constants.h`).
- Record: LCP byte, 10 packed k-mer bytes, mask byte, little-endian
  variable-width position then contig; **strand is the contig field's high bit**.
  Record buffers allocate `+1` byte for `msd_sort`'s sentinel — keep it.
- Forward packing covers `[j, j+40)` and stores `j`; reverse-complement covers
  `[j-28, j+12)` and stores `j+12`. `Get_Contig_Piece` writes a sentinel at
  `buffer[-1]`, so callers pass `buf+1`. Preserve syncmer selection and its
  tie behaviour (`mz == min4` falls through) exactly.
- "LCP" in this repo means **semantic bases-in-common in [0,40]**, not
  `msd_sort`'s raw byte 0 and not GIXshow's `lcp` column. `sort_records` alone
  does *not* recalculate it; `recalc_all_lcps` (called from `run_stage2`) does.
- Distributed sort = global prefix histogram → deterministic bucket ownership →
  `MPI_BYTE` `Alltoallv` → local radix sort → LCP recalculation → neighbour
  boundary check + seam fixup. Keep it deterministic and communication-free
  where it currently is.
- Contig-pair assignment: length-product weights, descending order, `(i,j)`
  tie-break, least-load heap with rank tie-break. `int8_t` storage caps world
  size at 127; the map is `ncontig²` bytes.
- Adapter structs mirror FastK/FastGA layouts and are **cast** to those
  interfaces — field order is load-bearing. Cumulative index uses `prefix <= p`
  semantics; clones share storage and only the wrapper is freed.
- A zero exit status does **not** mean PAF output was produced: `ALNtoPAF`
  failure is printed and then ignored.

## Routing

Read only what the task needs; these are not startup imports.

| Task | Guide |
|---|---|
| Targets, flags, stale-binary pitfalls | `docs/agent/build.md` |
| Phases, symbols, formats, MPI calls | `docs/agent/architecture.md` |
| salloc/srun, binding, paths, example jobs | `docs/agent/perlmutter.md` |
| `verify-phase*` prerequisites and commands | `docs/agent/validation.md` |
| Optimization proposal O1–O5, future tests | `docs/agent/research-plan.md` |
| Current objective, results, next action | `docs/agent/state.md` |

Nested notes: `src/AGENTS.md`, `scripts/AGENTS.md`, `lib/AGENTS.md`.

## Handoff convention

At the end of a work session, update `docs/agent/state.md` in place: current
objective, inspected commit, validation actually run (distinguish
source-inspected from executed), open unknowns, single next action. Keep it
short — it is a baton, not a changelog. Record newly found defects as
follow-up items there rather than fixing them mid-task.

Do not commit generated data, binaries, PAF/CSV outputs, or `tests/` fixtures.
