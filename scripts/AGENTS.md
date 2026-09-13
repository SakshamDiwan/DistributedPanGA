# scripts/ — verification and gold generation

Prerequisites, exact commands and current pass/fail status:
`docs/agent/validation.md`. Slurm details: `docs/agent/perlmutter.md`.

Three groups:

- **`verify_phase*.sh`** — regression checks, invoked via `make verify-phase1`,
  `verify-phase2`, `verify-phase3-stage1`, `verify-phase3-stage2`,
  `verify-phase3-stage2-lcp`. The last three need an allocation and honour
  `RANK_COUNTS` (default `1 2 4 8`). `verify_phase3_stage2_lcp.sh` starts with
  `make clean` — run it separately from any build you care about.
- **`make_gold.sh`, `build_lcp_oracle*.sh`, `normalize_gixshow*.c`,
  `diff_kmers.py`** — gold generation and oracles. `make_gold.sh` **deletes and
  rebuilds** `yeast7.1gdb`/`yeast7.gix` in place and needs `~/FASTGA`; treat
  running it as a deliberate act, never a side effect. Preserve existing gold
  artifacts and record their provenance.
- **`runtime_test_phase4b.sh`** — the only end-to-end driver. Not a `make`
  target, not parameterized by `RANK_COUNTS`.

Conventions to keep when editing:

- These scripts call `srun` themselves. Run them as plain commands inside an
  allocation; never `srun bash scripts/...`.
- `set -euo pipefail` throughout. Fixture checks come first and fail with a
  clear message, because a missing fixture is the common case here.
- Scripts key on dump columns by position: the default 4-column dump has the
  k-mer in `$4`, the `PGA_DUMP_LCP=1` 5-column dump has LCP in `$4` and k-mer
  in `$5`. Changing a dump format silently breaks the awk checks.
- Known non-portable assumptions (hardcoded `~/FASTGA`, another user's FastGA
  path, fixed `/tmp` filenames, CSV truncation) are documented in
  `docs/agent/perlmutter.md` and listed as follow-ups in
  `docs/agent/state.md`. Do not fix them as a drive-by.
