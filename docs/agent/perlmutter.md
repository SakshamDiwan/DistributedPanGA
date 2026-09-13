# Perlmutter operations

CPU account: **m4341**. `m4341_g` exists but is the GPU account and is out of
scope. Everything below was checked against the live system on 2026-09-10;
`[executed]` marks commands actually run, and no jobs were submitted.

## Verified account and hardware facts

`sacctmgr -n show assoc user=$USER account=m4341 format=account,qos` **[executed]**
lists, among others, `interactive`, `debug`, `jupyter`, `largemem`, `overrun`,
`preempt`, `cron`, plus the `gpu_*` set. So the `interactive` and `debug` QOS
are both available on m4341.

CPU node topology, from `scontrol show node` on an idle `regular_milan_ss11`
node **[executed]**:

```
Sockets=8  CoresPerSocket=16  ThreadsPerCore=2  CPUTot=256  RealMemory=515100
AvailableFeatures=cpu,milan,ss11
```

That is 128 physical cores / 256 logical CPUs and ~503 GiB usable per node.
Slurm reports 8 "sockets" because the node is in a 4-NUMA-per-socket
configuration; it is 2 physical AMD Milan sockets. `--constraint=cpu` selects
these nodes.

Re-check with `sinfo`/`scontrol` rather than trusting this file: NERSC changes
QOS names, limits and defaults, and the numbers above are a snapshot.

## Allocation

Interactive, one node, 30 minutes — the documented starting point:

```bash
salloc --account=m4341 --constraint=cpu --qos=interactive --nodes=1 --time=00:30:00
```

Roles:

- **`salloc`** — reserves nodes and drops you into a shell *on a login node*
  with the allocation attached. Use for iterating on validation. Commands you
  type run on the login node until you `srun` them.
- **`sbatch`** — queues a script; the script body runs on the first allocated
  node. Use for anything longer than an interactive session or that must
  survive a disconnect.
- **`srun`** — launches the parallel step *inside* an existing allocation. This
  is the only thing that creates MPI ranks. `mpirun` is not the supported
  launcher here; use `srun` even though some source comments say `mpirun`.

Status and control:

```bash
squeue --me                       # my queued/running jobs
squeue -j <jobid> -o "%.18i %.9P %.8T %.10M %.6D %R"
sacct -j <jobid> --format=JobID,JobName,State,Elapsed,MaxRSS,NNodes,NTasks
scontrol show job <jobid>         # full record incl. node list
scancel <jobid>                   # cancel
scancel -u $USER                  # cancel all of mine
```

For memory, `sacct ... MaxRSS` is per-step and needs the job to have finished;
state the measurement method whenever you report memory (see
`research-plan.md`).

## Shared versus node-local paths

| Path | Visibility | Use |
|---|---|---|
| `$HOME` (`/global/u2/...` = `/global/homes/...`) | shared, all nodes | source, small artifacts; **not** for heavy parallel I/O |
| `$SCRATCH` | shared, all nodes | GDB inputs and run outputs for real jobs |
| `$CFS/<project>` | shared | longer-lived datasets |
| `/tmp` | **node-local, per node** | where `pga-mpi-merge` puts its scratch |

This split matters because of two hardcoded behaviours:

1. `pga-mpi-merge` creates scratch with `mkdtemp("/tmp/pga-mpi-merge-r<rank>-XXXXXX")`.
   It is **node-local and `TMPDIR` is ignored** — setting `TMPDIR` will not
   relocate it. Ranks on different nodes cannot see each other's scratch. This
   is fine for the current design (scratch is never shared) but rules out any
   assumption that a rank can read a peer's `.1aln`.
2. `out.rank<N>.paf` and `pga-mpi-merge.timing.csv` are written to the process
   **current working directory** with relative paths. All ranks share that
   directory, so it must be on a shared filesystem for a multi-node run to
   collect every rank's PAF.

Rule: run the end-to-end binary from a **unique, shared, empty** working
directory, with an **absolute** GDB stem.

```bash
RUNDIR=$SCRATCH/pga/run-$(date +%Y%m%d-%H%M%S)-j${SLURM_JOB_ID:-manual}
mkdir -p "$RUNDIR" && cd "$RUNDIR"
```

A unique directory per run also keeps the appended timing CSV interpretable and
avoids mixing `out.rank*.paf` from different world sizes.

## PATH requirement

`pga-mpi-merge` converts its `.1aln` by shelling out to a **bare** `ALNtoPAF`,
so the FastGA tools must be on `PATH` in the environment `srun` propagates.
Verify before launching, not after:

```bash
command -v ALNtoPAF || { echo "ALNtoPAF not on PATH"; exit 1; }
export PATH=/path/to/FASTGA:$PATH
srun ... ./pga-mpi-merge "$GDB_STEM"
```

`ALNtoPAF` failure is only *printed*; `main` still returns 0. A successful exit
therefore does not mean PAF was produced — always check that
`out.rank*.paf` exist and are non-empty.

In this checkout `ALNtoPAF` is **not** on `PATH` and neither `~/FASTGA` nor
`/global/u2/s/sdiwan/FASTGA` exists **[executed]**, so no end-to-end run is
possible yet.

## Rank and thread placement

Two things consume cores per rank:

- `sort_records` passes **`nthreads = 4`** to `msd_sort` — hardcoded in
  `src/sort.c`, not configurable by env var.
- `NTHREADS = 1` in the merge/align drivers, so the FastGA phases are
  single-threaded per rank.

So budget **at least 4 physical cores per rank** or the sort threads will
oversubscribe. `OMP_NUM_THREADS` affects neither number; setting it is
misleading.

With 128 physical cores (256 logical) per node, and Slurm counting logical
CPUs in `-c`:

| Ranks/node | `-c` (logical) | Physical cores/rank | Note |
|---|---|---|---|
| 1 | 256 | 128 | 4 sort threads only |
| 2 | 128 | 64 | |
| 4 | 64 | 32 | |
| 8 | 32 | 16 | comfortable for 4 sort threads |
| 16 | 16 | 8 | |
| 32 | 8 | 4 | minimum for 4 sort threads |
| 64 | 4 | 2 | **oversubscribed** |

Use `--cpu-bind=cores` so the 4 sort threads land on distinct physical cores:

```bash
srun -N 1 -n 8 -c 32 --cpu-bind=cores ./pga-mpi-merge "$GDB_STEM"
```

Note `world_size > 127` is rejected outright by `contig_assignment_build`
(`int8_t` storage), so a single job cannot exceed 127 ranks today.

## Bounded example jobs

All are deliberately small. Do not scale these up without an explicit decision
about the experiment budget.

### Interactive smoke test, one node, 4 ranks

```bash
salloc --account=m4341 --constraint=cpu --qos=interactive --nodes=1 --time=00:30:00
# then, inside the allocation:
cd /global/u2/h/hsc53/DistributedPanGA
make CC=cc pga-mpi-merge
RUNDIR=$SCRATCH/pga/smoke-$SLURM_JOB_ID && mkdir -p "$RUNDIR" && cd "$RUNDIR"
export PATH=/path/to/FASTGA:$PATH
command -v ALNtoPAF || exit 1
srun -n 4 -c 64 --cpu-bind=cores \
     /global/u2/h/hsc53/DistributedPanGA/pga-mpi-merge /abs/path/to/yeast7
ls -l out.rank*.paf && cat pga-mpi-merge.timing.csv
```

### Batch, two nodes, 8 ranks

```bash
#!/bin/bash
#SBATCH --account=m4341
#SBATCH --constraint=cpu
#SBATCH --qos=debug
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=64
#SBATCH --time=00:20:00
#SBATCH --job-name=pga-mpi-merge
set -euo pipefail
REPO=/global/u2/h/hsc53/DistributedPanGA
export PATH=/path/to/FASTGA:$PATH
command -v ALNtoPAF >/dev/null || { echo "ALNtoPAF not on PATH" >&2; exit 1; }
RUNDIR=$SCRATCH/pga/run-$SLURM_JOB_ID
mkdir -p "$RUNDIR" && cd "$RUNDIR"
srun --cpu-bind=cores "$REPO/pga-mpi-merge" /abs/path/to/yeast7
ls -l out.rank*.paf
```

### Do not wrap a script that launches `srun` itself

`scripts/runtime_test_phase4b.sh` and the three MPI `verify_*` scripts call
`srun` internally. Run them as **plain shell commands inside** an allocation:

```bash
bash scripts/runtime_test_phase4b.sh        # correct
srun -n 8 bash scripts/runtime_test_phase4b.sh   # WRONG: 8 copies, each srun-ing again
```

## Known limitations of the existing scripts

Documented here rather than fixed; fixing them is separate work
(`state.md` follow-ups).

- `scripts/runtime_test_phase4b.sh` hardcodes `FASTGA_BIN_DIR=/global/u2/s/sdiwan/FASTGA`,
  a path that does not exist for this user. It also chooses its sweep from
  `$SLURM_JOB_NUM_NODES` (1, 2 or 4 only, erroring otherwise), deletes
  `out.rank*.paf` and truncates `pga-mpi-merge.timing.csv` at the start of every
  run, writes fixed `/tmp/gold.sorted.paf` and `/tmp/got.sorted.paf`, and runs
  from `$ROOT` — so its outputs land in the repository, not a run directory.
  `RANK_COUNTS` does **not** parameterize it.
- `scripts/make_gold.sh` assumes `$HOME/DistributedPanGA` and `~/FASTGA`.
- The `verify_*` scripts use `mktemp -d` under `/tmp`; on a multi-node
  allocation that directory is not shared, so the per-rank output files they
  concatenate must be produced by ranks on the same node as the script.
