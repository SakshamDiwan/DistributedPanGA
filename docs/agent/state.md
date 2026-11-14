# State / handoff

Last updated: 2026-10-12 (O1 closed).

## Current objective

**O1 is closed.** See `o1-final-report.md` (ACM PDF:
`docs/o1-report/o1-final-report.pdf`).
- The extraction speedup (O1a) is merged in PR #7.
- O1b, O1c, O1d and supermers were dropped, with evidence.

There is no active objective. The next candidate is alignment-stage
scheduling (proposal O3), pending a decision with the advisor
(`decision-gate-2.md`).

## Code status

| PR | content | state |
|---|---|---|
| #7 | week-1 verification, O1-0, O1a | merged |
| #8 | `ALNtoPAF` timeout, fatal on failure | merged |
| #9 | merge/alignment instrumentation, early memory frees | open, base `main` |
| #10 | exchanges > 2 GB per rank (bounded rounds) | open, stacked on #9 |
| #11 | `PGA_SCRATCH_DIR`, per-task cost diagnostics, decision gates | open, stacked on #10 |
| (this) | O1 final report, this file | stacked on #11 |

Merge in order: #9 → #10 → #11 → this.

## Environment

- Build and run under `module load cpu; module unload darshan` with `cc`
  (`build.md`).
- Run-time switches:
  - `PGA_PAIR_TASKS=1` writes the per-task alignment records
    (`pairs.rank*.tsv`); off by default.
  - `PGA_MEM_STAGE_RESET=0` keeps `/usr/bin/time`/`sacct` peak memory valid.
  - `PGA_EXCHANGE_ROUND_BYTES` sets the exchange round budget (default 1 GiB;
    the ranks use the minimum).
  - `PGA_SCRATCH_DIR` sets the scratch base (default `/tmp`, max 512
    characters).
  - `PGA_ALNTOPAF_TIMEOUT` sets the converter time limit (default 600 s).
- FastGA tools: `~/FASTGA` (`178afcb`).
- Datasets with FastGA gold: `$SCRATCH/pga-fixtures/`, with a `MANIFEST.txt`
  each. yeast7 is also symlinked from `tests/fixtures/`.

  | dataset | size | gold alignments |
  |---|---|---|
  | yeast7 | 83 Mbp | 43,126 |
  | scerevisiae8 | 96 Mbp | 50,046 |
  | primates4.chr6 | 679 Mbp | 1,547,485 |
  | primates14.chr6 | 2.64 Gbp | 4,577,556 |

  `$SCRATCH` is purged after ~8 weeks idle. Sources: yeast7 from
  `hypervolu.me/~erik/yeast`; the others from Zenodo 10.5281/zenodo.7937947.
- Job scripts: `jobs/` (copies of the scripts behind each result; see
  `jobs/README.md`). Logs stay in `~/pga-jobs/*/logs/`. The user submits all
  Slurm jobs.

## Validation executed (latest)

- **Correctness:** every completed run since PR #7 matches FastGA gold.
  scerevisiae8 at N=2 aborted on #9 alone (2 GB exchange limit); #10 fixes
  it. Runs covered:
  - yeast7, scerevisiae8: N = 1–8;
  - primates4: N = 4–16;
  - primates14: N = 16, and 32 on 2 nodes.
- **Test suites:** unit 44,231 checks plus ASan, and the MPI tier, on the
  #10/#11 code. On `followup/post-o1`: `make test` passes and the MPI tier
  passes at N = 1, 2, 3, 4, 8 with M-06 comparing all digests (job 60448488,
  `TMPDIR=$SCRATCH/tmp`).
- **Jobs:** 59712760, 59713118, 59715282, 59717508, 59717520. Results are in
  `pr1-results.md`, `pr2-results.md` and `decision-gate-{1,2}.md`.

## Open items (recorded, not fixed)

- **Extraction wait is timed as sort (found 2026-11-14):** each rank stops its
  extraction timer when it finishes (`pga-mpi-merge.c:324`, no barrier), then
  waits for the slowest extractor in the sort's first `MPI_Allreduce`
  (`dsort.c:384`). Stages are reduced separately with `MPI_MAX`, so the wait
  is counted twice: in `extract_s` (slowest rank) and in `sort_s` (fastest
  rank), and `total_s` overstates. Size = extraction spread across ranks:
  ~0 on yeast7 and primates4 at 4 ranks; ~2.6 s on primates4 at 8/16 ranks;
  4-9 s on primates14 (extraction splits whole contigs, so at 64 ranks 36 of
  64 ranks extract nothing). Report numbers are corrected in post-processing
  (`report/fig/make_figs.py`, `fix_sort_wait`: sort = max_r(extract_r +
  sort_r) - max_r(extract_r), from the `[o1-0]` lines); corrected totals agree
  with the stopwatch wall time to ~1 s. **Fix:** `MPI_Barrier` after
  `extract_to_records`, before `t_extract` is taken, so the wait counts as
  extraction. Runs made before the fix need the post-processing correction.
  Fix is on branch `timing-barrier` (26b80da, not merged; also logs
  `extract_wait_s` per rank and keeps `extract_s` as the rank's own time).
- **hsapiens90.chr6 at 64 ranks misses 1 of 1,375,248 alignments (found
  2026-11-14, open):** reproducible (jobs 61244457, 61277456); 32 ranks match
  exactly (61255434, 61277457); FastGA identical at -T16/-T32. Missing line:
  same-strand self-alignment of `HG03453#2#JAGYVV010000063.1` (GDB contig
  1072), a 6,419,539-6,423,436 vs b 6,419,557-6,423,438. `PGA_PAIR_TASKS`:
  pair (1072,1072,N) has the same 229 seeds and 6 hits at 32 and 64 ranks,
  but nlass 1 vs 0, so the difference is in alignment. Lead (not proven):
  `align_contigs` (identical to upstream FastGA 178afcb apart from debug
  code) handles a self-hit that straddles the main diagonal by setting only
  `path->abpos = path->aepos = 0`; `bbpos/bepos` keep the previous pair's
  values and feed `eant = aepos + bepos`, so later chains can be skipped.
  Trace build on branch `trace-align` (`PGA_TRACE_ANCHOR=1072`, job 61278174)
  shows the stale `b=87063320-87150831` from pair (1072,367). To prove: rerun
  the 32-rank trace (61278175 was cancelled), then test the one-line fix
  (also zero `bbpos/bepos`) at 32 and 64 ranks.
- **127-rank limit:** the contig-pair map is `int8_t`, and two `ncontig²`
  maps are live at once.
- **Large self-alignment tasks floor scaling:** primates14's largest is
  773 s. Splitting them while keeping FastGA-identical output is open.
- **Peak memory per stage:** with #9, stage 2 is the peak at 1 rank. O1c-style
  staging could lower it (optional).
- **Old defects still open:**
  - 128 MiB adapter index per rank;
  - dead `cleanup` label and `read_file_all`;
  - `O(ixlen)` `GoTo_Kmer_Index`;
  - script portability;
  - `Read_GDB` error check never fires (it returns 1, callers test `< 0`);
  - the CSV `total_s` is a sum of per-stage maxima.
- **M-06 can pass vacuously:** `tests/run_mpi.sh` skips missing digests and
  still prints "ALL MPI TESTS PASS". Under `salloc <cmd>` the script runs on
  the login node, so its `/tmp` work dir is invisible to the compute ranks.
  Fail on a missing digest and default the work dir to `$SCRATCH`.

## Next action

#9–#12 are in `main` (via #13–#15). Merge `followup/post-o1`. Then decide with the advisor whether to pursue O3,
a dynamic alignment-task scheduler with self-first ordering (modelled −33 to
−42% for the slowest rank at 16–32 ranks), or to stop here.
