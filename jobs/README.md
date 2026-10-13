# Slurm job scripts (Perlmutter, account m4341)

These are the exact scripts behind the results in `docs/agent/*.md`. They are
kept as a record, so they are not general-purpose tools.

## Path changes from the originals

- Absolute paths under the author's home directory were replaced by `$HOME`.
- Each `#SBATCH --output` now writes `%x.%j.log` to the submit directory.

## How they expect to be laid out

- The repository at `$HOME/DistributedPanGA`.
- FastGA binaries in `$HOME/FASTGA`.
- These scripts copied to `$HOME/pga-jobs/<dir>/`: several of them `source`
  a helper from there, e.g. `$HOME/pga-jobs/pr1/pr1_lib.sh`.
- Datasets in `$SCRATCH/pga-fixtures/`; see each `MANIFEST.txt`.

## Directories

| directory | purpose | results |
|---|---|---|
| `week1/` | week-1 verification and yeast7 baseline (A–D, C2) | `tests/RESULTS.md`, `docs/agent/baseline.md` |
| `o1/` | O1-0 instrumentation and O1a A/B runs | `docs/agent/o1-0-results.md`, `o1a-results.md` |
| `fix/` | ALNtoPAF timeout checks (R-05) | PR #8 |
| `data/` | FastGA gold for the Zenodo PGGB datasets (`make_gold_paf.sh`) | dataset `MANIFEST.txt` files |
| `pr1/` | instrumentation and memory validation | `docs/agent/pr1-results.md` |
| `pr2/` | large-count exchange validation | `docs/agent/pr2-results.md` |
| `step3/` | larger datasets, cost diagnostics, primates14 on 1 and 2 nodes | `docs/agent/decision-gate-{1,2}.md` |

## Safety notes

Every job checks its output against FastGA gold. Each timed run is wrapped in
an `srun -t` limit, so a hang is recorded instead of eating the allocation.
