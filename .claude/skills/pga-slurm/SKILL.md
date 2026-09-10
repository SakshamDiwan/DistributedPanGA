---
name: pga-slurm
description: Allocate, launch, monitor and cancel DistributedPanGA jobs on Perlmutter under CPU account m4341. Use when asked about salloc, sbatch, srun, rank and thread placement, CPU binding, job status, or where run outputs and scratch land.
---

Full reference, including verified node topology, path visibility rules and
bounded example jobs: `docs/agent/perlmutter.md`. Build first with
`docs/agent/build.md`.

## Workflow

1. Allocate. Interactive starting point:
   `salloc --account=m4341 --constraint=cpu --qos=interactive --nodes=1 --time=00:30:00`
   Use `sbatch` for anything longer or that must survive a disconnect. Confirm
   current QOS and limits with `sinfo`/`sacctmgr` rather than assuming.
2. Choose placement. Budget **at least 4 physical cores per rank**: the local
   sort hardcodes 4 threads (`src/sort.c`) and the merge uses `NTHREADS=1`.
   `OMP_NUM_THREADS` changes neither. CPU nodes have 128 physical / 256 logical
   cores, and Slurm's `-c` counts logical CPUs — e.g. 8 ranks/node → `-c 32`
   with `--cpu-bind=cores`. `world_size` above 127 is rejected by the code.
3. Prepare the working directory. `out.rank<N>.paf` and
   `pga-mpi-merge.timing.csv` are written to the **current working directory**
   by relative path, so use a unique, empty, **shared** directory under
   `$SCRATCH` and pass an **absolute** GDB stem. Internal scratch is a
   hardcoded per-rank `/tmp` directory and `TMPDIR` will not relocate it.
4. Verify `PATH` before launching: `command -v ALNtoPAF || exit 1`. The binary
   shells out to a bare `ALNtoPAF`, and a conversion failure is only printed —
   `main` still returns 0.
5. Launch with `srun` inside the allocation. Do **not** wrap a script that
   itself calls `srun` (the `verify_*` and `runtime_test_phase4b.sh` scripts
   do); run those as plain shell commands.
6. Monitor and clean up: `squeue --me`, `sacct -j <jobid> --format=...`,
   `scontrol show job <jobid>`, `scancel <jobid>`.
7. Confirm success by checking `out.rank*.paf` exist and are non-empty — not by
   the exit status.
