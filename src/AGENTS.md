# src/ — implementation

Full symbol-level map: `docs/agent/architecture.md`. Constraints summary:
root `AGENTS.md`.

Layering, innermost first: `constants.h` → `record`/`tables`/`pack`/`work_split`
→ `extract` → `sort` → `dsort` (MPI) → `kmer_adapter` → `contig_assignment` →
`fastga_pipeline` → drivers (`pga-*.c`).

- One driver per phase; `pga-mpi-merge.c` is the final one and defines the
  FastGA externs (`NTHREADS`, `FREQ`, `IBYTE`, `NCONTS`, `RankForPair`, …) that
  `fastga_pipeline.c` declares `extern`. Adding an extern means updating every
  driver that links the pipeline (`pga-merge.c` too).
- `fastga_pipeline.c` is lifted FastGA code, kept close to upstream on purpose.
  Match its brace and indentation style there rather than the surrounding code,
  and keep edits minimal so it stays diffable against FastGA.
- `kmer_adapter.h` structs are **cast** to `Kmer_Stream` / `Post_List`. Field
  order and widths are load-bearing — never reorder, insert, or change a type.
- `#include "constants.h"` defines `KMER 40`, which collides with the runtime
  `int KMER` in the merge drivers; `pga-mpi-merge.c` `#undef KMER` after its
  includes. Preserve that if you add includes.
- Buffers get `+1` byte for `msd_sort`'s sentinel; `Get_Contig_Piece` callers
  pass `buf+1` for its `buffer[-1]` sentinel. Both are easy to lose in a
  refactor and fail silently.
- `sort_records` hardcodes 4 threads; the merge drivers hardcode `NTHREADS=1`
  and `FREQ=10`. No env var overrides these.
- Errors in MPI code paths call `MPI_Abort` with a distinct code
  (2 = local sort check, 3 = boundary check). Keep codes distinct.

Headers are not Makefile prerequisites — after editing any `.h` here, rebuild
with `make -B` or the binary will be stale. See `docs/agent/build.md`.
