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
  #10/#11 code.
- **Jobs:** 59712760, 59713118, 59715282, 59717508, 59717520. Results are in
  `pr1-results.md`, `pr2-results.md` and `decision-gate-{1,2}.md`.

## Open items (recorded, not fixed)

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

## Next action

Merge #9–#11 and this PR. Then decide with the advisor whether to pursue O3,
a dynamic alignment-task scheduler with self-first ordering (modelled −33 to
−42% for the slowest rank at 16–32 ranks), or to stop here.
