---
name: pga-benchmark
description: Measure and record DistributedPanGA performance on Perlmutter with reproducible provenance. Use when asked to benchmark, time, profile, sweep rank counts, or interpret pga-mpi-merge.timing.csv stage timings.
---

Required record fields and the O1–O5 proposal context:
`docs/agent/research-plan.md`. Launching: `docs/agent/perlmutter.md`.
Correctness gate: `docs/agent/validation.md`.

## Workflow

1. Gate on correctness first. A performance number from an unvalidated build is
   not reportable — run the phase verifiers for the configuration you are about
   to measure.
2. Build cleanly and consistently: `make clean && make -B CC=cc all`. Record the
   expanded compile command from `make -Bn <target>` and confirm `-DLCPs`.
   Never mix `cc`- and `mpicc`-built binaries in one measurement.
3. Run from a unique, empty, shared working directory per configuration. The
   timing CSV is **appended**, so a shared directory mixes world sizes; the
   end-to-end script also truncates it and deletes `out.rank*.paf` on every
   invocation.
4. Measure **actual wall time separately**. The CSV's `total_s` is the sum of
   per-stage `MPI_MAX` values, not an end-to-end measurement, and it
   overestimates because different ranks can be the max in different stages.
   Note also that `merge` includes adapter setup and the seed exchange, and
   `alignments` includes the `ALNtoPAF` conversion.
5. Record every field listed in `research-plan.md`: input and reference hashes,
   source revision **and dirty state**, compiler/MPI versions and wrapper,
   flags, nodes/ranks/threads and binding, job ID, wall time separate from
   stage-max total, memory measurement method, repetitions with dispersion, and
   the correctness outcome.
6. Keep the baseline fixed — `FREQ=10`, `NTHREADS=1`, 4 sort threads — unless
   the experiment is explicitly about changing it, and say so if it is.
7. When comparing across world sizes, treat agreement as **consistency**
   evidence only. It is not an independent oracle; only FastGA-derived gold data
   is.
