# State / handoff

Last updated: 2026-11-23 (branch `r3-extract-sort`, not pushed; base `main`
5c1be75, PR #19 merged).

## Current objective

**R3 (extraction, distributed sort, record exchange): done on branch
`r3-extract-sort` (head after 93d70dd), report
`docs/agent/r3-extract-sort-report.md` (§0 summary, §1 checklist).** Not
pushed, no PR yet. Alignment code untouched; `pga-mpi-merge.c` changed only in
the extraction/sort block.

New defaults, each switchable back:
- extraction buffer freed after packing (`PGA_STAGE2_KEEP_LOCAL=1`);
- one `MPI_Alltoallv_c` for stage 2 (`PGA_EXCHANGE_LARGE_COUNT=0`);
- 1,024 sort parts, FastGA's layout (`PGA_SORT_PARTS=256`);
- heavy-part split (`PGA_SORT_HEAVY=0`).

Opt-in: `PGA_SORT_THREADS` (recommend `-c/2` in job scripts),
`PGA_SPLITTER=optimal` (not worth it). Barrier after extraction (26b80da) is
included. Effect: sort stage −22–29%, stage-2 per-rank memory −38% (human 64:
17.97 → 13.13 s, 14.07 → 8.66 GB).

Validation executed (2026-11-23, jobs 61880294–61883670, builds in
`~/pga-build/r3` = f57b432 and `~/pga-build/r3b` = 93d70dd; results in
`~/pga-r3-results/`):
- gold: 32/32 MATCH on yeast7, primates4, primates14, hsapiens90 at 1–64
  ranks;
- MPI tier PASS on all 3 builds;
- 357 `pga-dsort-bench` runs (8–256 ranks) give one index digest per dataset;
- `make test`, `make test-asan` pass.

**Single next action:** open the PR (formal, draft) and tell Saksham which
block of `pga-mpi-merge.c` changed; then optionally re-run the report sweep on
this branch for end-to-end numbers.

Open items found here:
- FastGA `msd_sort` threading defects (§6 of the report): threads never
  joined, `pthread_create` unchecked, `if` around `cond_wait`, static
  mutex/condvar re-initialised per call. Production is safe (one call per
  process); report upstream.
- ≥128 ranks: deeper prefix buckets and sub-contig extraction are deferred
  (§12.2).

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
  - `PGA_FASTGA_ORDER=0` turns off FastGA-order reproduction (needs
    `<stem>.gix`; without it the old order is used with a warning).
    `PGA_FASTGA_THREADS` (default 32) must match the thread count of the
    gold; `PGA_FASTGA_ORDER_LOG=1` logs every self pair's entry value.
- GitHub CLI: `gh` in `~/.local/bin`, logged in as the user; open PRs with
  `gh pr create --draft`.
- FastGA tools: `~/FASTGA` (`178afcb`).
- Datasets with FastGA gold: `$SCRATCH/pga-fixtures/`, with a `MANIFEST.txt`
  each. yeast7 is also symlinked from `tests/fixtures/`.

  | dataset | size | gold alignments |
  |---|---|---|
  | yeast7 | 83 Mbp | 43,126 |
  | scerevisiae8 | 96 Mbp | 50,046 |
  | primates4.chr6 | 679 Mbp | 1,547,485 |
  | primates14.chr6 | 2.64 Gbp | 4,577,556 |
  | hsapiens90.chr6 | 15.51 Gbp | 1,375,248 |

  `$SCRATCH` is purged after ~8 weeks idle. Sources: yeast7 from
  `hypervolu.me/~erik/yeast`; the others from Zenodo 10.5281/zenodo.7937947.
- Job scripts: `jobs/` (copies of the scripts behind each result; see
  `jobs/README.md`). Logs stay in `~/pga-jobs/*/logs/`. The user submits all
  Slurm jobs.

## Validation executed (latest)

- **2026-11-17, build 6032f52 (code unchanged since), jobs 61449300-61449335:**
  report sweep, 12 configurations x 3 repeats, all 36 gold=MATCH against the
  original FastGA gold; repeats within ~1%. Results in
  `~/pga-report-results-order/` (layout of `report/fig/make_figs.py`).
- **2026-11-16, jobs 61447693-61447707:** same 12 configurations once each,
  all MATCH; same-binary controls with `PGA_FASTGA_ORDER=0`: human 296.1 /
  186.0 s vs 301.5 / 188.6 s on (32 / 64 ranks; off loses the alignment at
  64), primates14 64 ranks 1,180.4 vs 1,188.3 s.
- **Values:** diagnostic FastGA (`~/pga-build/fastga-diag`, jobs
  61435904-61435906 + login-node yeast7) prints the entry `bepos` of every
  self pair; all 1,202 derived values match, and FastGA's thread ranges match
  ours.
- **Not executed on this branch:** `make test` / MPI test tier; no unit
  tests exist for `fastga_order.c`.
- Older (O1-era) validation: jobs 59712760, 59713118, 59715282, 59717508,
  59717520; `pr1-results.md`, `pr2-results.md`, `decision-gate-{1,2}.md`.

## Open items (recorded, not fixed)

- **Extraction wait is timed as sort (found 2026-11-14):** fixed by the
  barrier (26b80da), now on `r3-extract-sort` (0714a54), not yet on `main`.
  Runs made without it need `report/fig/make_figs.py`'s `fix_sort_wait`.
- **hsapiens90.chr6 at 64 ranks missed 1 of 1,375,248 alignments (found
  2026-11-14; fixed on branch `match-fastga-order`, 2026-11-16, not merged):**
  - Cause (traced, jobs 61278174 / 61295918): `align_contigs`, identical to
    FastGA 178afcb `FastGA.c:3257`, skips the aligner for a same-strand self
    hit that straddles the main diagonal and zeroes only `abpos/aepos`; the
    loop then reads `path->bepos` left by the previous real aligner call
    (`eant = aepos + bepos`) and jumps if it exceeds the current position.
    FastGA's pair order fixes that value; ours depended on rank ownership
    (653 at 32 ranks, 87,150,831 at 64 for contig 1072).
  - Fix chosen: reproduce FastGA's value (`src/fastga_order.c`,
    `Self_Order` hooks in `src/fastga_pipeline.c`). FastGA order = contig
    permutation in `<stem>.gix` (length-sorted), parts by length, thread
    ranges by seed count (`PGA_FASTGA_THREADS`, default 32). Each
    earlier-partner pair's final bepos is sent to the owner of the self pair;
    self pairs wait only for the answers they need. `PGA_FASTGA_ORDER=0`
    restores the old path. Needs `<stem>.gix`; adds an N*N byte map.
  - Executed (jobs 61447693-61447707, one run each, final build 6032f52):
    gold=MATCH on yeast7 (1,2,4,8 ranks), primates4 (4,8,16), primates14
    (16,32,64), human (32,64). Same binary, feature on vs off: human 301.5 vs
    296.1 s (32), 188.6 vs 186.0 s (64, off loses the alignment), primates14
    64 ranks 1,188.3 vs 1,180.4 s. Values: all 1,202 derived entry values
    equal those printed by a diagnostic FastGA (`~/pga-build/fastga-diag`,
    jobs 61435904-61435906 and a login-node yeast7 run).
  - Not derivable: where FastGA has no earlier call in a thread range it
    reads an unset field (11 self pairs over the four datasets). Observed
    32764-32767 for the first contig of the run, 0 elsewhere; we use 32767
    (`PGA_FASTGA_UNSET_BEPOS`) and 0.
  - Unit tests (2026-11-18): `make test-fastga-order` (F-01..F-06,
    `tests/unit/test_fastga_order.c`) on the MPI-free core
    `src/fastga_order_core.c`; ranges against a verbatim copy of FastGA's
    loops, all rules against a sequential model of FastGA incl. self pairs
    depending on self pairs (up to 11 rounds). Nine mutations all caught.
    Full `make test` passes. The MPI messaging/scheduling is still covered
    only end to end.
  - Load balancing (O3), replay of recorded per-pair times
    (`scripts/simulate_pair_schedule.py`, records of jobs 59717520 and
    61277456): with the current length-product assignment the order costs
    0%. With a cost-based assignment it costs 28% on primates14 at 16 ranks
    (contig 0 is short but its self pair is the most expensive, 773 s, and its
    0.008 s prerequisite runs late on another rank); ordering equally near
    earlier partners by fewest seeds first removes that (0% everywhere
    tested). A dynamic queue costs up to 69% unless prerequisites of
    expensive self pairs get priority, then 0%. Human: 0% for every static
    variant, up to 24% for the dynamic queue without priority. Splitting a
    self pair along diagonals: only the chunk containing the main diagonal
    reads the leftover (code reading, not tested).
  - Open: the self-depends-on-self round
    never occurred in a dataset run; the off path matched gold on yeast7 (1 process)
    and on human 32 / primates14 64 only; three repeats of the headline
    timings on this build; report to FastGA's authors. The alternative
    one-line fix (zero `bbpos/bepos`) is on local branch `fix-self-straddle`:
    it changes up to 0.053% of FastGA's own alignments, so it was not taken.
  - **Report sweep on this build (2026-11-17, jobs 61449300-61449335):** the
    report's 12 configurations, three repeats each, all 36 gold=MATCH;
    repeats within about 1%. Filed in `~/pga-report-results-order` (layout of
    `report/fig/make_figs.py`); tables, figures and every number in the
    report text now come from it. `make_figs.py` reads hsapiens90 from
    results.tsv, leaves `table_fastga.tex` alone when there is no
    `fastga.tsv` (FastGA rows are still Saksham's and the hand-added human
    ones), and subtracts each rank's wait for other ranks from its
    alignment time.
  - **Unexplained, for Saksham:** his earlier table rows for the three small
    datasets were 4-20% slower in total than this sweep (yeast7 1 rank
    79.97 vs 66.61 s; per-rank search 41.91 vs 23.08 s). Not a code effect:
    the old code in our job setup already gave the new times (yeast7 8 ranks
    9.40 s, primates4 16 ranks 70.5 s on 2026-11-14; yeast7 1-rank search
    22.6 s on 2026-10-12). His FastGA times do reproduce here within 2%.
    So his pipeline runs differ in setup (build flags, binding, darshan?);
    his `jobs/report` scripts are not in the repo. Ours:
    `jobs/human/07_report_sweep.sbatch`.
  - User decides when to merge. Offline notes (not in git):
    `notes/local/missing-alignment.pdf`, `replicating-fastga.pdf`.
- **Report (2026-11-14, branch `report-human-merge`, not pushed):** human rows
  are now in the shared tables (`report/fig/table_{pipeline,fastga,perrank}.tex`,
  added by hand; Saksham's generator is not in the repo and would overwrite
  them) and both figures are redrawn by `report/fig/make_report_figs.py` from
  those table files. The report leads with the 32-rank human result (matches
  FastGA). Still to do: pipeline figure, extraction speedup numbers, 3 repeats
  and 16-ranks-per-node placement for human, FastGA at 64 threads on human.
- **Aligned bases (gold PAFs, PAF column 11 summed):** yeast7 0.55 Gbp,
  primates4 4.67, primates14 137.9, hsapiens90 642.7. Human has 4.7x the
  aligned bases of primates14 in 0.30x the alignments, at about half the
  summed alignment time; why is untested.
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
