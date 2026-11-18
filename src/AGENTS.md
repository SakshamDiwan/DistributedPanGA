# src/ — implementation

Full symbol-level map: `docs/agent/architecture.md`. Constraints summary:
root `AGENTS.md`.

Layering, innermost first: `constants.h` → `record`/`tables`/`pack`/`work_split`
→ `extract` → `sort` → `dsort` (MPI) → `kmer_adapter` → `contig_assignment` →
`fastga_pipeline` → drivers (`pga-*.c`).

- One driver per phase; `pga-mpi-merge.c` is the final one and defines the
  FastGA externs (`NTHREADS`, `FREQ`, `A_PAYLOAD_BYTES`, `NCONTS`, `RankForPair`, …) that
  `fastga_pipeline.c` declares `extern`. Adding an extern means defining it
  in `pga-mpi-merge.c`.
- `fastga_pipeline.c` is lifted FastGA code. Keep its **behaviour** identical to
  upstream; its **names** are being made descriptive on purpose (e.g. `IBYTE` →
  `A_PAYLOAD_BYTES`). Pure renames must compile to identical machine code
  (compare `.text`/`.rodata`/`.data` of the two builds) and never touch log
  strings. Match its brace and indentation style there rather than the
  surrounding code.
- `fastga_order.c` + the `Self_Order` hooks in `fastga_pipeline.c` make each
  same-strand self pair start from the `path->bepos` FastGA's own pair order
  would have left (FastGA reads it after skipping a diagonal self hit). This
  is what makes the output independent of rank count; do not reorder the
  same-strand passes or drop the hooks without rerunning the gold comparison
  at several rank counts. `PGA_FASTGA_ORDER=0` gives the old single pass.
- `kmer_adapter.h` structs are used directly by `fastga_pipeline.c` (the
  `First_Kmer_Entry`-style macros call the `adapter_*` functions, no casts), so
  the compiler checks every field use. Never reintroduce a cast to
  `Kmer_Stream` / `Post_List`: that would make field order load-bearing again.
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
