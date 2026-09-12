# Test results

Last run: 2026-09-12, Perlmutter login node (gcc 14.3.0, cray-mpich).
Command: `make test` — no allocation, no fixtures, no FastGA tools. ~25 s.

**42,471 checks · 0 failures · 5 expected failures · exit 0**

## Unit suites

| Binary | Covers | Checks | Result |
|---|---|---|---|
| `test_record_split` | E-12, E-13, R-03 sizing | 112 | pass, 3 xfail |
| `test_extract` | E-01…E-11, R-03 extraction | 7,519 | pass, 1 xfail |
| `test_extract_window` | E-11 (`SCAN_MAX=64`) | 19 | pass |
| `test_lcp_sort` | S-01…S-05b, R-04 | 4,028 | pass, 1 xfail |
| `test_dsort` | S-06…S-09 | 26,668 | pass |
| `test_assign` | A-01…A-05 | 76 | pass |
| `test_adapter` | K-01…K-05 | 4,034 | pass |
| `test_exchange_plan` | R-02 | 15 | pass |

`make test-asan` (ASan + UBSan) also clean, and is the only run that includes
**S-05b**, which extracts 1.5 Mbp to force a production buffer growth.

## Known defects

| ID | What | Status |
|---|---|---|
| R-02 | Cumulative send displacement could overflow `int` | **fixed** — guard added in `src/seed_exchange.c`, markers removed |
| R-03 | `post_bytes` one byte short at `maxctg = 256^k`; a 256-base contig's RC record at position 256 decodes as **0** | open — 4 xfail assertions |
| R-04 | `recalc_all_lcps` leaves a singleton's LCP byte stale (observed 37, contract says 0) | open — 1 xfail assertion |
| R-01 | LCP seam computed against an empty rank's sentinel instead of the true predecessor | open — test written, needs an allocation to run |

## Not yet run

The MPI tier (**M-01…M-06**, **R-01**) is written and compiles, but multi-rank
runs need a Slurm allocation:

```bash
salloc --account=m4341 --constraint=cpu --qos=interactive --nodes=1 --time=00:30:00
make test-mpi                     # RANK_COUNTS="1 2 3 4 8"
```

Single-rank smoke test already passes on the login node (`uniform` and `skew`
modes, world=1). R-01 needs exactly 4 ranks.

**P6 end-to-end (Z-01, Z-02) is blocked**: no `tests/fixtures/`, no `~/FASTGA`,
no `ALNtoPAF` on `PATH`.

## Notes

- Tests build with `$(CFLAGS)`, so `-DLCPs` matches the `LCPS=1` default.
  Without it `recalc_all_lcps` silently takes its 0/1 branch.
- MPI binaries need `MPICH_GPU_SUPPORT_ENABLED=0` in the default Perlmutter
  environment, which loads `craype-accel-nvidia80`; otherwise cray-mpich aborts
  with *"GPU_SUPPORT_ENABLED is requested, but GTL library is not linked"*.
  `tests/run_mpi.sh` sets it.
- Two fixture assumptions were confirmed empirically rather than assumed:
  `select[768]=0, select[769]=1, select[770]=2` (R-01's empty middle rank), and
  that removing the `+1` from the production allocators makes S-05b fail under
  ASan.

See `tests/CATALOG.md` for what each test asserts.
