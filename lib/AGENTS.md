# lib/ — vendored FastGA/FastK compatibility boundary

`GDB.c`, `gene_core.c`, `libfastk.c`, `ONElib.c`, `MSDsort.c`, `align.c`,
`alncode.c`, `ANO.c`, `RSDsort.c` are **upstream FastGA/FastK sources**, kept
here so every target compiles them directly (there is no library archive).

Treat these as read-only. Our code is written to match them, not the reverse:

- `src/kmer_adapter.h` structs mirror `Kmer_Stream` (`libfastk.h`) and are cast
  to it; `src/fastga_pipeline.h`'s `Post_List` mirrors FastGA's. Any field
  reordering here breaks those casts silently.
- `msd_sort` (`MSDsort.c:404`) writes a sentinel at `array[nelem*rsize]`, which
  is why record buffers allocate `+1`, and its internal byte 0 encoding is
  **not** semantic LCP — see `docs/agent/architecture.md`. `rmsd_sort`
  (`RSDsort.c:292`) is used by the seed sort in `fastga_pipeline.c`.
- `Get_Contig_Piece` (`GDB.h`) writes a sentinel at `buffer[-1]`, so callers
  pass `buf+1`.
- The build is warning-noisy from these files (`align.c`, `ANO.c`); that is
  pre-existing. Only new warnings from `src/` are findings.

If a change genuinely requires touching `lib/`, record which upstream file and
revision diverged in `docs/agent/state.md` — otherwise the next FastGA update
will silently revert or conflict with it.
