# R3 — extraction, distributed sort and record exchange

Branch `r3-extract-sort` (base `main` 5c1be75, PR #19 merged). Started 2026-11-23.
Scope: extraction, prefix-bucket partitioning, record exchange, local sort,
stage-2 memory. Alignment scheduling, pair ownership (`contig_assignment.c`),
the FastGA-order mechanism and `fastga_pipeline.c` are untouched.

Status legend: **DONE** (measured or implemented and tested here),
**MEASURED / NO CHANGE NEEDED**, **DEFERRED**, **BLOCKED**.

**Outcome in one paragraph.** The stage-2 (extract → sort) path is now
22–29% faster and needs 38% less memory per rank, with byte-identical output.
Three new defaults do it:
1. the extraction buffer is freed after packing;
2. the local sort uses FastGA's own 1,024-part layout and splits heavy
   (poly-A/T) parts;
3. the exchange is one `MPI_Alltoallv_c`.

Human 64 ranks: sort stage 17.97 → 13.13 s, per-rank peak 14.07 → 8.66 GB.
The summed per-rank peak, compared with FastGA, falls from 49.6× to 29.0×
(human 32 ranks: 44.8× → 27.6×). End-to-end gold matched FastGA in all 32
runs (yeast7, primates4, primates14, hsapiens90; 1–64 ranks; defaults and
variants). The sorted index (k-mers and LCPs) is identical across 357 sort
benchmarks at 8–256 ranks. The five-base bucket imbalance is real but not the
bottleneck below ~128 ranks; the local sort was, and inside it one thread
sorting the AAAAA bucket.

All numbers below come from scripts in the repo, run on data that is on disk:

| source | what | how to regenerate |
|---|---|---|
| `~/pga-report-results-order/` | the 36-run report sweep (2026-11-17, jobs 61449289–61449326, build 6032f52, all gold=MATCH), per-rank `[o1-0]`, `[pr1]`, `[pr1-mem]` lines | `scripts/r3_stage2_breakdown.py` |
| `pga-prefix-hist` output (scratch) | exact 12-base prefix histogram + per-contig record counts of hsapiens90, primates14, primates4, yeast7, through the production extraction scan | `make pga-prefix-hist`; `scripts/r3_bucket_analysis.py` |
| `report/fig/r3_*.pdf` | figures | `scripts/r3_figs.py --hist <dir>` (needs `module load python`) |
| `~/pga-r3-results/` | new benchmark and gold runs (jobs 61880294–61883670, 2026-11-23): round 1 `<ds>/` + `gold.tsv` (opt-in variants, build f57b432), round 2 `defaults2/` (1,024 parts + large-count default, cc9c5bc), round 3 `bench2/` + `defaults3/` (heavy-part split, eec599b/93d70dd) | `jobs/r3/submit_r3.sh`, `scripts/r3_bench_summary.py <dir>` |

Validation of the offline method: the record totals from `pga-prefix-hist`
equal the logged totals (human 11,439,544,067), and the Python port of
`compute_destination_map` reproduces the **logged per-rank slice sizes
exactly** for every configuration in the sweep (human 32/64, primates14
16/32/64, primates4 4/8/16); the port of `compute_work_split` reproduces the
logged extracted-records-per-rank exactly. So the histograms and simulations
below are the production distribution, not an approximation.

Caveat on the sweep logs: they predate the `timing-barrier` fix, so each
rank's `hist_s` contains its wait for the slowest extractor. All sort-stage
numbers below remove it (`sort_wall = max_r(extract_r + sort_r) - max_r(extract_r)`,
the same correction `report/fig/make_figs.py` applies). The fix itself is now
on this branch (cherry-picked 26b80da).

---

## 1. R3 checklist

| item | status | one-line result |
|---|---|---|
| Timing: extraction wait counted as sort | **DONE** | `timing-barrier` (26b80da) cherry-picked; new `[r3-s2]` line splits every substage into local work vs. collective wait |
| P0-A bucket histogram, per-rank loads | **DONE** | buckets are very skewed (max/mean 8.97, CV 0.92) but per-rank slices are balanced to 1.05 (32 ranks) / 1.18 (64 ranks); the 64-rank excess comes from the greedy splitter, not from indivisible buckets; the 5-base scheme hits a hard floor at ≥128 ranks |
| P0-B sort scaling 32→64 | **DONE** | before: 22.1 → 16.8 s (66%) on the *same 4 nodes*. Limiting component: the local sort (56%), above all one thread sorting the AAAAA bucket; then the node-bandwidth-bound exchange (21–25%). After the changes: 17.95 → 13.13 s, and 128 ranks on 8 nodes 8.34 s |
| P0-C memory | **DONE** | stage-2 peak = extraction + receive + (send or sort scratch) + 2 GiB staging, model within 0.3 GB of all 36 sweep runs; the 45–50× is human-only (sum over ranks of in-memory index copies vs. FastGA's out-of-core GIXmake) |
| Memory: free extraction buffer after packing; drop staging | **DONE, default** | per-rank peak human 32: 22.88 → 14.21 GB; 64: 14.07 → 8.66 GB; aggregate human 32: 675 → 416 GB (44.8× → 27.6× FastGA) |
| P1-A partitioning alternatives | **DONE** (analysis + opt-in optimal splitter, measured) | any split that keeps whole 12-base groups is exact; the optimal contiguous splitter helps ≤ 0.4 s at 64 ranks and hurts at 256, so it stays opt-in; deeper (≤12-base) prefixes are the fix for ≥128 ranks (DEFERRED) |
| P1-B radix sort | **DONE** (documented; 1,024 parts + heavy-part split default) | FastGA's in-place MSD radix sort; FastGA's GIXmake gives it 1,024 parts and we gave it 256; heavy parts now split; msd_sort 8.35 → 4.42 s at human 64 |
| P1-C extraction imbalance | **MEASURED / NO CHANGE NEEDED (≤64 ranks)** | whole-contig split: human max/mean 1.12 / 1.35 at 32 / 64; primates14 36 of 64 ranks idle; costs ≤ 1.5 s of a 186–1,180 s run |
| P2-A overlap | **MEASURED / NO CHANGE NEEDED** | the record exchange is ≤ 4.6 s, ≤ 2.2% of end-to-end on human and < 0.1% on primates14; not worth an asynchronous redesign |
| P2-B chunked exchange | **DONE** (`MPI_Alltoallv_c` default) | rounds were forced by 32-bit counts; staging cost 2 GiB/rank and two copies; `MPI_Alltoallv_c` is faster in every configuration (human 32: 4.63 → 3.27 s; primates14 16: 1.65 → 0.60 s) |
| P1-D 128+ ranks | **DONE** (measured at 128/256 ranks with `pga-dsort-bench`) | sort side runs at 128 and 256 ranks; 256 is not faster than 128 on 8 nodes (8.79 vs 8.34 s): five-base buckets and the remaining heavy groups bound it; the 127 cap is the `int8_t` pair map (alignment side, not touched) |
| N1 lemma | **DONE** (§9) | checked against `new_self_merge_thread` line by line; the exact condition is "whole 12-base groups per rank", weaker than "whole 5-base buckets" |
| End-to-end gold on this branch | **DONE** | 32/32 gold=MATCH: 15 (round 1), 7 (round 2), 10 (round 3, final defaults + physical-core threads) |
| MPI test tier (M-07, large-count) | **DONE** | PASS at 1, 2, 3, 4, 8 ranks on each of the three builds; M-06 made all 4 comparisons |
| FastGA `msd_sort` defects found | **RECORDED** (not fixed: `lib/` is upstream) | never joins its threads, ignores `pthread_create` failure, re-initialises its static mutex/condvar per call; safe for one call per process, which is how production uses it (§6) |

---

## 2. Current implementation (code map)

| step | code | notes |
|---|---|---|
| work split | `src/work_split.c:compute_work_split` | whole contigs; rank r gets contigs whose cumulative bases cross r·total/P |
| extraction | `src/extract.c:scan_contig`, `extract_scan` (new, exposes the scan with a callback), `extract_to_records` | single thread per rank; record buffer starts at 2^20 records and doubles |
| record | `src/record.h` | 1 LCP + 10 k-mer + 1 mask + position (LE) + contig/strand (LE); 16 B yeast7, 17 B primates, 18 B human |
| bucket | `src/dsort.h:record_bucket` | first 5 bases = 10 bits = 1,024 buckets |
| histogram | `dsort.c:bucket_local_counts`, `MPI_Allreduce` | 8 KiB per rank |
| ownership | `dsort.c:compute_destination_map` (FastGA `distribute()` lines 678–697); new opt-in `compute_destination_map_optimal` (`PGA_SPLITTER=optimal`) | contiguous bucket ranges, computed redundantly on every rank |
| pack | `dsort.c:organize_send_buffer` | two passes, out-of-place into a send buffer grouped by destination |
| exchange | `chunked_exchange.c:chunked_alltoallv` | default (R3): one `MPI_Alltoallv_c` when both sides are contiguous buffers; otherwise (seed exchange, or `PGA_EXCHANGE_LARGE_COUNT=0`) one `MPI_Alltoallv` if totals ≤ budget and ≤ INT_MAX, else rounds through staging |
| local sort | `sort.c:sort_records` → `lib/MSDsort.c:msd_sort` | out-of-place regroup into 1,024 parts (R3 default; was 256), heavy parts regrouped into 256 sub-parts, then one in-place MSD radix sort call, one thread per part; `PGA_SORT_THREADS`, `PGA_SORT_PARTS`, `PGA_SORT_HEAVY` |
| LCP | `sort.c:recalc_all_lcps`; seam fix in `run_stage2` | semantic LCP in bases |
| checks | `verify_per_rank_sorted`, `verify_cross_rank_boundaries` | `MPI_Abort` 2 / 3 |
| driver | `pga-mpi-merge.c` lines ~330–375 | now calls `run_stage2_take` and prints `[r3-s2]` |

New instrumentation (all off the fast path except ~6 reads of
`/proc/self/status` per run): `[r3-s2] rank r:` with `count_s allreduce_s
map_s pack_s sizes_s xfer_s free_s rebucket_s msd_s lcp_s check_s
boundary_s`, bytes sent/self/received, rounds, staging, sort settings,
parts, splitter, and `rss{0..5}_kb hwm{0..5}_kb` at substage ends.
`PGA_BUCKET_HIST=<file>` dumps the global histogram and bucket→rank map.

---

## 3. P0-A — prefix-bucket imbalance

Figures: `report/fig/r3_bucket_hist.pdf`, `r3_rank_loads.pdf`,
`r3_balance_vs_ranks.pdf`. Tables: `<tag>.buckets.tsv`, `<tag>.P<P>.rankloads.tsv`,
`<tag>.schemes.tsv` from `r3_bucket_analysis.py`.

### 3.1 The 1,024-bucket histogram is very uneven

| | hsapiens90.chr6 | primates14.chr6 | primates4.chr6 |
|---|---|---|---|
| records | 11,439,544,067 | 1,938,068,671 | 495,810,447 |
| non-empty buckets | 1,024 | 1,024 | 1,024 |
| max bucket | 100,227,207 (AAAAA, 0.876%) | 16,542,766 (AAAAA, 0.854%) | 4,378,293 (AAAAA, 0.883%) |
| min bucket | 1,180 (CCGTA) | 220 | 50 |
| mean / median | 11.17 M / 8.58 M | 1.89 M / 1.46 M | 0.484 M / 0.374 M |
| max / mean | 8.97 | 8.74 | 9.04 |
| coefficient of variation | 0.918 | 0.914 | 0.919 |
| next largest | TTTTT 0.83%, ATTTT 0.44%, TTTTG 0.39%, AAAAG 0.39% | same set | same set |
| largest 12-base group | TTTTTTTTTTTT 20.1 M (0.175%) | 0.170% | AAAAAAAAAAAA 0.176% |

Two sequence effects explain the shape: poly-A/T runs (AAAAA and TTTTT alone
hold 1.7% of all records, ~9× the mean) and CpG depletion (every one of the
smallest buckets contains `CG`). The shape is nearly identical across the
three primate datasets, so it is a property of primate DNA and of the
syncmer selection, not of the assembly.

### 3.2 Per-rank loads under the production splitter

The splitter assigns contiguous bucket ranges at the cumulative quantiles,
so many small buckets average out and only the largest buckets matter.

hsapiens90.chr6 (`prod-5base-greedy` = production):

| ranks | max/mean | CV | empty ranks | ranks < 0.5× mean | largest bucket / ideal rank load |
|---|---|---|---|---|---|
| 32 | 1.054 | 0.032 | 0 | 0 | 0.28 |
| 64 | 1.177 | 0.059 | 0 | 0 | 0.56 |
| 128 | 1.217 | 0.098 | 0 | 0 | **1.12** |
| 256 | 2.243 | 0.238 | 2 | 3 | 2.24 |
| 512 | 4.486 | 0.442 | 16 | 42 | 4.49 |

primates14 and primates4 give the same picture (64 ranks: 1.159 / 1.155;
128 ranks: 1.31 / 1.24).

Measured slice imbalance in the sweep equals these (human: max slice 6.78 GB
at 32 ranks, 3.79 GB at 64).

### 3.3 Where the imbalance comes from

| human | production | optimal contiguous, same 1,024 buckets | LPT, 1,024 buckets, non-contiguous | optimal contiguous, 6-base | 8-base | 12-base |
|---|---|---|---|---|---|---|
| 32 ranks | 1.054 | 1.027 | 1.000 | 1.008 | 1.001 | 1.000 |
| 64 ranks | 1.177 | 1.052 | 1.001 | 1.016 | 1.002 | 1.000 |
| 128 ranks | 1.217 | 1.122 | 1.122 | 1.034 | 1.003 | 1.000 |
| 256 ranks | 2.243 | 2.243 | 2.243 | 1.274 | 1.009 | 1.001 |

- **Up to 64 ranks the imbalance is the splitter's**, not the buckets': the
  largest bucket is ≤ 0.56 of a rank's ideal load, and the optimal contiguous
  split of the same buckets gets to 1.05. The FastGA-style greedy rule
  (`compute_destination_map`) rounds each cut to the nearer bucket edge and
  never revisits it, so errors accumulate.
- **From ~128 ranks the buckets themselves are the floor**: AAAAA alone is
  1.12× (128) and 2.24× (256) of a rank's ideal load, and no assignment of
  5-base buckets can do better. Beyond 1,024 ranks some ranks must be empty.
- Is it meaningful *today*? At 64 ranks the slowest rank sorts and receives
  17.7% more than the mean. Local sort and exchange time follow the largest
  slice, so the optimal splitter would cut them by up to ~11% (predicted
  ~1–1.5 s of a 16.8 s sort stage, < 1% of the 186 s end-to-end). Real but
  small. Measured (§7.2): the optimal splitter saves 0.43 s at human 64 and
  costs 1.5 s at 256 ranks, so it stays opt-in. The local sort's internal
  balance (one thread per part, §6) turned out to matter far more than the
  slice balance across ranks.

---

## 4. P0-B — why the sort scales weakly from 32 to 64 ranks

Figure: `report/fig/r3_stage2_breakdown.pdf`. Medians of 3 repeats from the
sweep. Columns are the slowest rank per substage; "exchange" is measured on
the slowest packer (the others wait in the size `Alltoall` for it), "verify"
on the slowest sorter. The five components sum to the measured stage wall
time within 0.4–1.8%, so the decomposition is complete.

| config | sort wall (s) | hist | pack | exchange | local sort + LCP | verify | exchange share | local-sort share | moved | max slice |
|---|---|---|---|---|---|---|---|---|---|---|
| human 32 r / 4 n | 22.14 | 0.37 | 3.16 | 4.64 | 12.35 | 1.42 | 21% | 56% | 205.9 GB | 6.78 GB |
| human 64 r / 4 n | 16.78 | 0.23 | 2.27 | 4.17 | 9.41 | 0.72 | 25% | 56% | 205.9 GB | 3.79 GB |
| primates14 16 r / 1 n | 9.38 | 0.20 | 1.81 | 1.72 | 4.84 | 0.73 | 18% | 52% | 32.9 GB | 2.11 GB |
| primates14 32 r / 2 n | 6.49 | 0.17 | 1.45 | 1.16 | 3.20 | 0.41 | 18% | 50% | 32.9 GB | 1.08 GB |
| primates14 64 r / 4 n | 4.76 | 0.13 | 1.20 | 0.92 | 2.26 | 0.19 | 19% | 48% | 32.9 GB | 0.60 GB |

(The sweep's raw `sort_s` was 24.7 / 19.4 s for human; the difference is the
extraction wait. The earlier draft's 22.65 / 18.16 s are close to the raw
numbers of an older build; use the corrected 22.1 / 16.8 s.)

Strong-scaling efficiency of the sort stage: human 32→64 **66%**;
primates14 16→32 72%, 32→64 68%, 16→64 49%.

Per component, human 32→64 (ideal 2×): verify 1.98×, hist 1.6×, pack 1.39×,
local sort 1.31×, **exchange 1.11×**.

Why:

1. **Both human runs use the same 4 nodes.** Going from 32 to 64 ranks adds no
   network injection bandwidth and no memory bandwidth. The exchange moves the
   same 205.9 GB, ~3/4 of it between nodes (uniform destinations): about
   154 GB / 4 nodes / 4.2 s ≈ 9 GB/s per node, including two staging copies
   per byte (§8). It cannot shrink with more ranks per node; it shrinks with
   more nodes (primates14 1→2→4 nodes: 1.72 → 1.16 → 0.92 s).
2. **Local sort follows the largest slice**, which shrinks 1.79×, not 2×,
   because the splitter's max/mean goes from 1.054 to 1.177 (§3).
3. **Sort throughput per rank drops** (median 44.7 → 36.2 M records/s):
   - each rank's slice covers fewer first-byte values, and `msd_sort`
     parallelises only across those parts (§6). Modelled efficiency of 4
     threads over the actual part sizes: mean 0.79 → 0.67, worst rank
     0.41 → 0.25 (a rank whose slice is a single first byte sorts on one
     thread);
   - twice as many sorting threads per node (8 → 16 ranks × 4 threads)
     share the same memory bandwidth.
   primates14 runs at a constant 16 ranks/node and still loses throughput
   per rank (31.6 → 29.6 → 25.9 M rec/s) as parts per rank fall
   (16.7 → 8.7 → 4.7), which points at (3a). The bench settles it (§7.4):
   more threads alone change nothing (human 64: 8.36 → 8.34 s), 1,024 parts
   cut `msd_sort` by 35% and splitting the heavy part by another 18%, so
   (3a) dominated; (3b) shows up only as a smaller gain from extra threads.
4. Collective waits are small: `hist` + `verify` are ≤ 8% of the stage.
   Computation (pack + local sort + LCP + checks) is ~70–75% of the stage
   and communication 18–25%.

**The limiting component is the local sort** (~56%), then the
node-bandwidth-bound exchange (~21–25%). The optimisations in §7 target
exactly these two.

Requested-per-configuration fields not in the table (per-rank means, bytes
sent/received per rank, peak records per rank) are in
`r3_stage2_breakdown.py --tsv` output (`*_mean`, `*_max`, `recv_max_GB`,
`slice_max`). For the new runs `[r3-s2]` gives `sent_bytes`, `self_bytes`,
`recv_bytes` per rank directly.

---

## 5. P0-C — memory

### 5.1 What is live in stage 2 (per rank)

With L = extracted bytes on the rank, R = received slice bytes and S = exchange
staging (2 GiB when the chunked path runs, 0 otherwise), the old code holds:

| moment | live buffers |
|---|---|
| exchange | extraction buffer L + send buffer L + receive buffer R + staging S |
| local sort | extraction buffer L + receive buffer R + regroup scratch R |

Model: peak = L + R + max(L, R) + S. Against all 36 sweep runs the model is
within 0.3 GB of the measured stage-2 VmHWM wherever stage 2 is the peak
(human 32: 23.26 predicted vs 22.87 GB measured; human 64: 14.30 vs 14.06;
primates14 64: 6.84 vs 6.73). The 2 GiB staging was the unexplained 1.8–2 GB
residual of a model without S; yeast7 (single `Alltoallv`, S = 0) has none.

Other observations:
- **Extraction capacity vs. utilisation:** the record buffer doubles, so it is
  50–100% used (human 32 ranks: ~357 M records in a 537 M-record buffer).
  Untouched pages are never resident, and the measured extraction peak equals
  L to within 0.01–0.15 GB, so the unused capacity costs address space, not RSS.
- **`msd_sort` allocates nothing**: it permutes in place (cycle-leader with a
  20-deep stack). The regroup scratch in `sort_records` is ours.
- **LCP** costs nothing extra (byte 0 of each record).
- **Adapter** (next stage): recv R plus the adapter's own copy (13/18 of R)
  plus a 128 MiB prefix index (2^24 × 8 B) per rank: 12.4 GB at human 32.
- `PGA_MEM_STAGE_RESET` (default on) resets VmHWM at every stage boundary, so
  the per-stage peaks are separate high-water marks, **not additive**. The
  per-rank peak is the maximum over stages (`overall_peak_kb`); sacct/`time -v`
  then show only the last stage. Aggregates below sum each rank's
  `overall_peak_kb`, which is an upper bound on simultaneous use. It is tight
  for stage 2 because stage 2 is synchronous.

### 5.2 Per-stage peaks (max over ranks, GB; sweep, old code)

| config | extract | stage 2 | adapter | rest (merge/redist/align) | per-rank peak | sum over ranks |
|---|---|---|---|---|---|---|
| human 32 | 7.19 | **22.87** | 12.42 | < 3.2 | 22.87 | 675 |
| human 64 | 4.41 | **14.06** | 7.06 | < 3.1 | 14.06 | 747 |
| primates14 16 | 3.76 | **11.52** | 3.90 | — | 11.52 | 141 |
| primates14 64 | 2.46 | **6.73** | 1.34 | — | 6.73 | 238 |
| primates4 16 | 1.35 | **4.76** | 1.12 | — | 4.76 | 53 |
| yeast7 1 | 1.02 | 3.03 | 1.97 | align 5.32 | 5.32 | 5.3 |

Stage 2 is the per-rank peak in every multi-rank configuration.

### 5.3 The 45–50× comparison with FastGA

FastGA peaks (`/usr/bin/time`, same nodes): human `FastGA -T32` 15.06 GB
(job 61244455); `GIXmake -T16` 4.70 GB (job 61244453); primates14 `FastGA
-T32` 55.6 GB (job 61435905); primates4 14.0 GB (job 61435904).

| dataset, ranks | sum of rank peaks | FastGA peak | ratio |
|---|---|---|---|
| human 32 | 675 GB | 15.06 GB | **44.8×** |
| human 64 | 747 GB | 15.06 GB | **49.6×** |
| primates14 16 / 32 / 64 | 141 / 188 / 238 GB | 55.6 GB | 2.5× / 3.4× / 4.3× |
| primates4 4 / 8 / 16 | 33 / 43 / 53 GB | 14.0 GB | 2.4× / 3.1× / 3.8× |

So the reviewer's 45–50× is the human dataset. Its causes, largest first:
1. **Design:** the distributed index lives in RAM. Human has 11.44 G records ×
   18 B = 206 GB of records, which the old stage 2 holds about three times
   over across the ranks (extraction, send/scratch, receive). FastGA's
   GIXmake sorts out of core: it writes sorted index parts to disk and peaks at
   4.7 GB. FastGA's alignment then needs 15 GB, which is small on human
   (1.4 M alignments) but 55.6 GB on primates14, hence the much smaller ratio
   there.
2. **Per-rank fixed costs** grow with rank count: 2 GiB staging + 128 MiB
   adapter index + MPI buffers, so 64 ranks cost ~72 GB more than 32.
3. The per-stage resets make "sum of per-rank peaks" an upper bound, but stage
   2 is synchronous, so for stage 2 the bound is close to the true
   simultaneous use.

Per rank, DistributedPanGA at 32–64 ranks (14–23 GB) is about the same as
FastGA's single process (15 GB) on human. The paper should report both views
(see §11).

### 5.4 Change made, and projections

`run_stage2_take` frees the extraction buffer immediately after packing (it is
dead once copied into the send buffer). New model: peak = max(L + R, 2R) + S.
With `PGA_EXCHANGE_LARGE_COUNT=1` S = 0 as well (§8).

| config | old measured | model: take | model: take + large-count | next stage's peak |
|---|---|---|---|---|
| human 32 | 22.87 | 16.1 | 13.9 | adapter 12.4 |
| human 64 | 14.06 | 10.0 | 7.8 | adapter 7.1 |
| primates14 16 | 11.52 | 11.0 (adapter-bound 10.95) | — | — |
| primates14 64 | 6.73 | 5.2 | 5.2 | — |
| primates4 16 | 4.76 | 3.9 | 3.9 | — |

Aggregate for human (sum of rank peaks): 675 → ~486 GB at 32 ranks (44.8× →
~32×), 747 → ~567 GB at 64 (49.6× → ~38×); with large-count a further
~2 GB × ranks.

**Measured** (bench `s2_hwm_kb` and gold `[pr1-mem]` agree to 0.01 GB):

| config | pre-R3 | take only | take + large-count (default) | sum over ranks of max(extract, stage 2), pre → default | × FastGA peak |
|---|---|---|---|---|---|
| human 32 | 22.88 | 15.85 | **14.21** | 675 → 416 GB | 44.8 → 27.6 |
| human 64 | 14.07 | 9.82 | **8.66** | 748 → 436 GB | 49.6 → 29.0 |
| human 128 | 9.46 | 6.72 | **5.64** | 880 → 467 GB | — |
| human 256 | 6.82 | 4.92 | **4.62** | 994 → 540 GB | — |
| primates14 16 | 11.52 | 7.79 | **7.48** | 129 → 76 GB | 2.3 → 1.4 |
| primates14 64 | 6.74 | 4.76 | **4.77** | 158 → 113 GB | 2.8 → 2.0 |

Per-rank columns are the stage-2 peak in GB. The sums cover extraction and
stage 2 only; for human, stage 2 is the whole-run peak, so they equal the
end-to-end sums. For primates the alignment stage of FastGA and of ours
dominates both sides, so compare whole-run peaks there. The same-binary
control (`keep`, gold job 61880973) gave 22.87 GB at human 32. Stage 2 is still the per-rank peak on human; on primates14 64 the overall
peak is now another stage (5.25 GB end-to-end). The summed peak still grows
with rank count (fixed per-rank memory: adapter index, MPI, GDB), so ranks
should not be added for memory reasons alone.

Projection to larger inputs (records per base = 0.737 on all three primate
sets; record size 18 B while ≤ 2^15 contigs and contigs ≤ 4 Gbp): one copy of
the index is 13.3 B/bp in aggregate. Old stage 2 ≈ 40 B/bp + 2 GiB/rank;
new ≈ 26.5 B/bp + 2 GiB/rank (26.5 B/bp with large-count). Example: a
100 Gbp pangenome at 128 ranks needs ~20.7 GB per rank in stage 2 (new) vs
~31 GB (old), i.e. ~330 GB vs ~500 GB per node at 16 ranks/node; the old code
would not fit 512 GB nodes once staging and the adapter are added.

Further reductions (not done; each is a real change, not a buffer-lifetime
fix):
- Pack in place (permute the extraction buffer by destination, American-flag
  style) and regroup in place before `msd_sort`: removes one R-sized copy;
  peak becomes L + R in both phases. Gain on top of the change above: none at
  the exchange (send + receive are still both live), so not worthwhile alone.
- Stream the exchange in rounds that free send data and sort completed parts
  as they arrive: peak → ~R + one round. This is the only route to ≈1× the
  index in aggregate (~206 GB human). It interacts with overlap (§8) and is
  the right follow-up only if memory becomes the binding constraint.
- Out-of-core index as in GIXmake: removes the RAM requirement, at the cost of
  the I/O the in-memory design avoids.

---

## 6. P1-B — the radix sort

- **Origin:** `lib/MSDsort.c` is FastGA's `MSDsort.c`, vendored unchanged
  (`lib/AGENTS.md`). FastGA's GIXmake uses it to sort each index part.
- **Why MSD:** `radix_sort` starts at the most significant key byte (`digit`)
  and distributes records into 256 sub-buckets by that byte, then recurses
  into each sub-bucket on the next byte. Leading bytes shared by a whole
  bucket are skipped without a pass. Buckets of ≤ 15 records go to a shell
  sort. The distribution is **in place**: per byte value it keeps
  `off[x]`/`end[x]` and moves records along permutation cycles through a
  20-deep stack of slots (American-flag sort), so it needs no second array.
  It is not stable.
- **Keys:** `ksize = 1 + KBYTES = 11`, starting at `digit = 1`, so the keys are
  record bytes 1–10, the packed 40-mer, MSB first. Byte 0 is the LCP slot, and
  the mask, position and contig bytes are payload. Records with equal k-mers
  come out in an unspecified order (§7.3).
- **Pre-bucketing and `part[]`:** `msd_sort(array, nelem, rsize, ksize, part,
  beg, end, nthreads)` expects the array already grouped into contiguous parts
  whose concatenation, once each part is sorted, is sorted. `part[x]` is the
  **byte length of part x**. Before R3, `sort_records` built that layout by
  counting records per first k-mer byte (256 parts), copying them out of place
  into part order, and copying back; it now does the same with the 1,024
  five-base buckets.
- **Threading:** `msd_sort` keeps a stack of free thread ids and, for each
  non-empty part in index order, takes one (waiting on a condition variable
  when all are busy) and runs `pthread_create(sort_thread)` on that part. A
  part is sorted entirely by one thread. Parallelism is therefore min(threads,
  non-empty parts), and a dominant part serialises the sort.
- **Why LCPs are recomputed:** `msd_sort` writes byte 0 as it goes, using one
  formula inside a part (digit·4 + `LCP_Table` of the differing byte) and
  another at part seams (only when parts share bits 8–9). Neither is "bases in
  common with the predecessor" everywhere. `recalc_all_lcps` overwrites every
  byte 0 with `simple_kmer_lcp_bases`, and `run_stage2` fixes each rank's
  first record against the last record of the nearest non-empty preceding rank.
- **Four threads:** hard-coded until now. In the report placement (`-c = 256 /
  ranks-per-node`, `--cpu-bind=cores`) a rank owns 16 physical cores at 32
  ranks on 4 nodes and 8 at 64, so the sort uses 25% / 50% of the allocated
  cores. And a rank's slice spans only ~256/P first-byte values (8.8 parts at
  32 ranks, 4.8 at 64, 2.7 at 128), so even 4 threads are often starved:

  | human, 4 threads | 32 ranks | 64 ranks | 128 ranks |
  |---|---|---|---|
  | 256 parts: mean / worst-rank efficiency | 0.79 / 0.41 | 0.67 / 0.25 | 0.45 / 0.25 |
  | 1,024 parts: mean / worst-rank efficiency | 0.90 / 0.63 | 0.83 / 0.42 | 0.73 / 0.25 |

  (List-scheduling model over the real part sizes, work ∝ records;
  `<tag>.msdsched.tsv` also has 2/8/16 threads.)
- **FastGA itself uses 1,024 parts.** GIXmake calls `msd_sort(sarray, next,
  swide, KBYTES, panel, kbeg, kend, NTHREADS)` with `panel[p]` over the
  `NUM_BUCK = 1024` five-base buckets of its part (`notes/GIXmake_annotated.c:1449–1474`).
  Our 256-part regroup is the deviation, and the comment in `sort_records`
  claiming it "matches FastGA's internal layout" was wrong.
- **Now the default: 1,024 parts** (`PGA_SORT_PARTS=256` restores the old
  layout). This regroups by the dsort bucket, FastGA's layout. It is safe with
  the unmodified `msd_sort`: its seam code reads `LCP_Table[x ^ n]` only when
  parts x and n share bits 8–9, so x ^ n < 256; those byte-0 values are
  overwritten anyway. More than 1,024 part indices would read past
  `LCP_Table`.
- **Now the default: heavy-part split** (`PGA_SORT_HEAVY=0` disables it). A
  part holding more than 1/threads of the rank's records (and ≥ 65,536) is
  regrouped stably into 256 sub-parts on the next 8 key bits. All light parts
  and sub-parts are then merged, in key order, into ≤ 1,001 contiguous
  segments for **one** `msd_sort` call. Segment boundaries are key boundaries,
  so the result is sorted and the (k-mer, LCP) sequence is unchanged (SO-06,
  and the bench `kseq` digests). Measured: on human the AAAAA bucket
  (100,227,207 records) was the slowest rank's only large part at 64, 128 and
  256 ranks and took 5.2–5.5 s on one thread. The split cuts the slowest
  rank's `msd_sort` to 4.42 / 3.19 / 3.37 s. What remains is AAAAAAAAA, 33% of
  the bucket, as predicted from the histogram.
- `PGA_SORT_THREADS=n` uses n threads (default 4). Alone it does nothing,
  because parts, not threads, were the limit. Combined with 1,024 parts and
  the split it gains 2–6% (human 64: 13.13 → 12.81 s). It is not the default,
  because the right value is the rank's physical cores (`-c / 2` under
  `--cpu-bind=cores`) and the code cannot know the placement. Job scripts
  should set it; gold matched with it (`phys`).
- **Defects in FastGA's `msd_sort`** found while testing (`lib/MSDsort.c`,
  upstream, not changed):
  - it creates one thread per part and never joins or detaches them, so each
    exited thread keeps its stack mapping;
  - it ignores `pthread_create`'s return value;
  - it waits with `if (Tavail <= 0) pthread_cond_wait(...)`, not `while`;
  - it re-initialises its static mutex and condition variable on every call,
    while the previous call's last worker may still be inside
    `pthread_cond_signal`.

  Together they make a process hang once it has created about 32 k threads
  (`vm.max_map_count` = 65,530 mappings; observed in a unit test at exactly
  65,530). Back-to-back calls are also unsafe. Production calls `msd_sort`
  once per process with ≤ 1,024 threads, so it is unaffected. The heavy-part
  path keeps to one call for this reason, and the sort unit tests run each
  test in a forked child. Worth reporting to FastGA's authors.

---

## 7. Changes implemented and their measured effect

All on `r3-extract-sort`. Each default can be switched back by environment
variable, so the pre-R3 path is still available in the same binary.

| change | default | switch back | correctness argument | tests |
|---|---|---|---|---|
| `MPI_Barrier` after extraction (26b80da) | on | — | timing only | — |
| `[r3-s2]` per-rank substage line; `PGA_BUCKET_HIST` dump | on / off | — | read-only | M-07 checks the byte counters |
| `run_stage2_take`: extraction buffer freed after packing | **on** | `PGA_STAGE2_KEEP_LOCAL=1` | the buffer is never read after `organize_send_buffer`; same bytes sent | M-07 (slices byte-identical to `run_stage2`, input emptied) |
| exchange by one `MPI_Alltoallv_c` (contiguous buffers only) | **on** | `PGA_EXCHANGE_LARGE_COUNT=0` | same counts and displacements, 64-bit; seed exchange (file callbacks) unchanged | test_chunked_mpi large-count cases (MPI tier) |
| 1,024 sort parts (FastGA's layout) | **on** | `PGA_SORT_PARTS=256` | parts are contiguous key ranges in order | SO-01..05 (+ASan), M-07 |
| heavy-part split, one `msd_sort` call | **on** | `PGA_SORT_HEAVY=0` | sub-parts and segments are contiguous key ranges in order | SO-06 (+ mutation), bench digests |
| `PGA_SORT_THREADS=n` | 4 | — | threads sort disjoint parts | SO-01, SO-04/05 |
| `PGA_SPLITTER=optimal` | off | — | same invariants as the greedy map | S-10, M-07 |
| `pga-dsort-bench`, `pga-prefix-hist`, `r3_*.py` | — | — | tools | digests in every bench run |

### 7.1 Benchmark design

`jobs/r3/r3_dsort_bench.sbatch`: `pga-dsort-bench` = extraction, barrier,
stage 2, barrier, digests. Variants are interleaved inside one allocation,
3 repeats each, in the report placement (`-c = 256 / ranks-per-node`,
`--cpu-bind=cores`). Tables report medians. The min–max spread over the 3
repeats was ≤ 1% of the median in 101 of 118 cells and ≤ 2.4% in all
(`r3_bench_summary.py --tsv` has min and max). Sort wall = time between the
two barriers.

### 7.2 Round 1 — each option on its own (build f57b432; old defaults; jobs 61880928–61880939)

Sort-stage wall time (s):

| config | old | take | 1,024 parts | phys. threads | 1,024 + threads | large-count | 256 MiB rounds | 4 GiB rounds | optimal split | all |
|---|---|---|---|---|---|---|---|---|---|---|
| human 32 r / 4 n | 23.37 | 23.32 | 20.61 | 23.18 | 19.03 | 21.89 | 23.16 | 23.46 | 23.18 | 17.71 |
| human 64 r / 4 n | 17.59 | 17.58 | 14.61 | 17.56 | 14.72 | 16.60 | 17.40 | 17.69 | 17.15 | 13.84 |
| human 128 r / 4 n | 14.21 | 14.16 | 14.37 | 14.19 | 14.35 | 13.20 | 14.13 | 14.21 | 14.20 | 13.40 |
| human 128 r / 8 n | 11.08 | 11.13 | 11.07 | 11.10 | 11.07 | 10.45 | 11.10 | 11.08 | 11.09 | 10.45 |
| human 256 r / 8 n | 11.55 | 11.56 | 11.62 | 11.54 | 11.63 | 10.91 | 11.48 | 11.56 | 13.04 | 12.49 |
| primates14 16 r / 1 n | 9.81 | 9.81 | 9.11 | 9.26 | 8.39 | 8.69 | 9.74 | 9.94 | 9.82 | 7.37 |
| primates14 32 r / 2 n | 6.76 | 6.76 | 5.84 | 6.73 | 5.56 | 6.19 | 6.67 | 6.95 | 6.60 | 5.05 |
| primates14 64 r / 4 n | 4.99 | 4.99 | 4.19 | 5.00 | 4.17 | 4.64 | 5.00 | 5.11 | 5.00 | 3.90 |
| yeast7 8 r / 1 n | 0.40 | 0.41 | 0.39 | 0.33 | 0.32 | 0.40 | 0.41 | 0.41 | 0.40 | 0.31 |

Slowest rank's `msd_sort` and exchange transfer (s):

| config | msd: 256 parts | 1,024 | 1,024 + threads | exchange: rounds (1 GiB) | `Alltoallv_c` |
|---|---|---|---|---|---|
| human 32 | 9.61 | 6.92 | 5.50 | 4.63 | 3.27 |
| human 64 | 8.36 | 5.43 | 5.41 | 3.87 | 2.94 |
| human 128 / 8 n | 5.24 | 5.22 | 5.21 | 2.39 | 1.76 |
| human 256 / 8 n | 5.47 | 5.47 | 5.48 | 2.72 | 2.04 |
| primates14 16 | 3.59 | 2.81 | 2.08 | 1.65 | 0.60 |
| primates14 64 | 2.10 | 1.30 | 1.29 | 0.89 | 0.55 |

Findings:
- **take** costs no time and saves 31–38% per-rank memory.
- **1,024 parts** is the largest single gain at ≤ 64 ranks.
- **Extra threads** help only together with more parts.
- **`Alltoallv_c`** wins everywhere, both by dropping the staging copies and
  through the library's own algorithm.
- **Round budget:** 256 MiB is marginally faster and lighter than 1 GiB; 4 GiB
  is slower and heavier.
- **The optimal splitter** gains ≤ 0.4 s at 64 ranks and loses 1.5 s at 256.
  There the AAAAA bucket fixes the max, so min-max filling pushes many ranks
  up to that max and leaves 129 empty.
- From 64 ranks up, `msd_sort` stays at 5.2–5.5 s whatever the setting. The
  slowest rank is the one holding the AAAAA bucket as a single part, which led
  to the heavy-part split.

### 7.3 Round 3 — final defaults (build eec599b/93d70dd; jobs 61883641–61883670)

`pre` = every R3 change off; `noheavy` = final defaults without the heavy-part
split; `phys` = final defaults + `PGA_SORT_THREADS` = physical cores per rank.

| config | pre | noheavy | **default** | phys | default vs pre | msd: pre → default | per-rank peak: pre → default |
|---|---|---|---|---|---|---|---|
| human 32 r / 4 n | 22.99 | 19.04 | **17.95** | 16.94 | −22% | 9.62 → 5.90 | 22.88 → 14.21 GB |
| human 64 r / 4 n | 17.97 | 14.02 | **13.13** | 12.81 | −27% | 8.35 → 4.42 | 14.07 → 8.66 GB |
| human 128 r / 8 n | 11.09 | 10.48 | **8.34** | 8.01 | −25% | 5.25 → 3.19 | 9.46 → 5.64 GB |
| human 256 r / 8 n | 11.54 | 10.95 | **8.79** | 8.78 | −24% | 5.46 → 3.37 | 6.82 → 4.62 GB |
| primates14 16 r / 1 n | 9.81 | 8.02 | **8.03** | 7.04 | −18% | 3.59 → 2.81 | 11.52 → 7.48 GB |
| primates14 64 r / 4 n | 4.86 | 3.73 | **3.43** | 3.26 | −29% | 2.10 → 0.98 | 6.74 → 4.77 GB |
| yeast7 8 r / 1 n | 0.41 | 0.39 | **0.39** | 0.32 | −5% | 0.15 → 0.13 | 0.42 → 0.29 GB |

Sort-stage strong scaling, human, default: 32 → 64 ranks (same 4 nodes)
68% (pre: 64%); 64 ranks / 4 nodes → 128 ranks / 8 nodes 79%; 128 → 256 ranks
on 8 nodes none (8.34 → 8.79 s), the five-base floor of §3.

### 7.4 End-to-end correctness and effect

32 gold runs, all `gold=MATCH` (sorted byte-for-byte PAF against the original
FastGA gold, duplicates kept):

| round | build | variants | configurations | result |
|---|---|---|---|---|
| 1 | f57b432 | `default` (take), `fast` (1,024 parts + phys. threads + large-count + optimal split), `keep` | yeast7 1/8, primates4 4/16, primates14 64, human 32/64 (+ `keep` human 32) | 15/15 MATCH |
| 2 | cc9c5bc | new defaults (1,024 parts + large-count) | same 7 | 7/7 MATCH |
| 3 | 93d70dd | final defaults (+ heavy split); `phys` | same 7; `phys` on primates14 64, human 32/64 | 10/10 MATCH |

Multi-rank MPI tier: PASS on all three builds (1, 2, 3, 4, 8 ranks, incl. M-07
and the large-count cases; M-06 made 4 comparisons). Sort-index digests: one
(`kseq`, `multi`) pair per dataset over 357 bench runs (human 198, primates14
114, yeast7 45), across every variant and 8–256 ranks.

End-to-end sort stage (gold jobs, one run each, barrier build):

| config | pre-R3 sweep (corrected) | final default | phys | per-rank peak (whole run) |
|---|---|---|---|---|
| human 32 | 22.14 | 18.26 | 17.18 | 22.87 → 14.21 GB |
| human 64 | 16.78 | 13.16 | 13.29 | 14.06 → 8.66 GB |
| primates14 64 | 4.76 | 3.58 | 3.33 | 6.73 → 5.25 GB |

The end-to-end wall time is dominated by alignment, whose run-to-run spread
on these single runs (±1–2%) is as large as the sort-stage saving (4–5 s on
human, 1 s on primates14). Claim the saving for the sort stage and the memory
reduction, not for end-to-end time.

## 8. P2 — overlap and the chunked exchange

### 8.1 Overlap (P2-A)

Upper bound on what overlap could hide = the exchange time (assuming perfect
overlap and no slowdown of the overlapped computation):

| config | exchange | sort stage | end-to-end | max saving / end-to-end |
|---|---|---|---|---|
| human 32 | 4.64 s | 22.1 s | 301 s | 1.5% |
| human 64 | 4.17 s | 16.8 s | 187 s | 2.2% |
| primates14 64 | 0.92 s | 4.8 s | 1,176 s | 0.08% |

Feasibility of each candidate:
- *Extraction vs. exchange:* ownership needs the global histogram, which needs
  all of extraction. Overlap would need an ownership map from a sample, and a
  fallback for misestimates. Not worth 1–2%.
- *Packing vs. communication:* the bounded-round structure already allows it
  (pack round k+1 while round k is in flight). Packing is 2.3 s at human 64.
  The gain is bounded by min(pack, exchange).
- *Receiving vs. sorting:* parts arrive mixed from every source, so a part is
  complete only when all sources have delivered it. Sorting completed parts
  early needs sources to send in bucket order (they could, since the send
  buffer can be bucket-ordered at no cost) and rounds aligned to part
  boundaries. Ordering and ownership would be unchanged, and the sorted output
  identical.

**Recommendation:** do not build an asynchronous exchange for speed. If memory
becomes binding, the streaming receive-and-sort design (§5.4) is worth doing
for memory, and overlap would come with it.

### 8.2 Chunked exchange (P2-B)

- **Why rounds:** `MPI_Alltoallv` takes `int` counts and displacements in
  units of the datatype (`MPI_BYTE`). Human slices are 3.8–6.8 GB per rank,
  above 2^31 − 1, so a single call is impossible. The code moves the same
  bytes in rounds of at most `budget / P` per peer, which fits `int`
  displacements, through staging buffers. Human 32: 8 rounds; human 64: 5;
  primates14 64: 3; yeast7: 1 (direct path, no staging).
- **Memory:** staging is bounded by budget (default 1 GiB) per direction: 2
  GiB per rank whenever chunking happens, on top of the send and receive
  buffers. Bounded, but not free; it is the 1.8–2 GB residual of the memory
  model.
- **Above 32-bit limits:** displacements within a round are bounded by budget
  (≤ INT_MAX enforced by `xchg_chunk_cap`); totals are `int64`. Correct at any
  size (R-06 tests mixed budgets down to one record per peer).
- **`MPI_Alltoallv_c` (now the default):** cray-mpich 9.1.0 is MPI 4.1 and provides it.
  `MPI_Count` counts and `MPI_Aint` displacements remove the rounds, the 2 GiB
  staging and the two staging copies (each byte is copied into send staging and
  out of receive staging). It does **not** remove the send or receive buffers.
  Used only when both sides are contiguous (stage 2), so the file-streamed
  seed exchange is unchanged; `PGA_EXCHANGE_LARGE_COUNT=0` restores the rounds.
  Measured (§7.2): faster in every configuration (human 32: 4.63 → 3.27 s;
  human 128 on 8 nodes: 2.39 → 1.76 s; primates14 16: 1.65 → 0.60 s), and
  −2 GiB per rank.
- **Budget** (with rounds): 256 MiB is 0.1–0.2 s faster than 1 GiB and saves
  1.6 GB per rank; 4 GiB is 0.1–0.4 s slower and costs 2 GB more. Rounds are
  not free; the copies through staging are what they cost.
- **Extra copies avoidable without MPI-4?** Only by sending directly from the
  send buffer with per-round displacement offsets, which `int` displacements
  forbid once the buffer exceeds 2 GiB. A derived datatype
  (`MPI_Type_contiguous` of a record) would raise the limit 18× without
  staging, but the counts must then be in records. That is a smaller change
  than staging and works on MPI-3 too. Not done; large-count is simpler where
  available.

---

## 9. N1 — communication-free seed generation (lemma and proof)

### 9.1 Facts used, checked against the code

1. *Adapter groups.* `build_kmer_adapter` (`src/kmer_adapter.c:38–93`) splits a
   rank's sorted slice into **prefix groups** of equal 24-bit prefix = the
   first 12 bases (`prefix_bytes = 3`), in slice order, and copies each
   record's LCP byte verbatim.
2. *Merge reads one group at a time.* `new_self_merge_thread`
   (`src/fastga_pipeline.c:244–405`, FastGA's merge) copies exactly one prefix
   group into `cache`, writes a marker LCP of 11 after its last entry, and
   processes only that cache until it is exhausted. All seed pairs it emits
   are (suf1, p) with both entries in `[low, hgh)` inside the cache, i.e. in
   the same 12-base group.
3. *Frequency test.* For each entry the merge finds the current LCP run
   `[low, hgh)` at length `plen ≥ 12`, and emits seeds only if `hgh < low +
   FREQ` entries (`kfreq = FREQ·entry_bytes`). `low` and `hgh` are positions
   in the cache, so the run is counted within the group, and the group is all
   of the run's occurrences on this rank.
4. *Use of LCP bytes.* Inside the group the merge reads LCP bytes of entries
   2..n (each ≥ 12) and the marker (11). For the first entry it reads only
   `rcur[LCP]`, as the lower bound of `for (i = rcur[LCP]; i <= plen; i++)
   vlcp[i] = rcur`. The only `vlcp` index ever read below 12 is 11 (a
   single-entry group, where `plen = 11`). So every value of the first entry's
   LCP in [0, 11] produces identical behaviour. Its true value is < 12 by
   definition of a group start. `vlcp[]` persists between groups, so it also
   matters that no stale entry is read. At refill `vlcp[12]` and
   `vlcp[first-LCP..plen]` are written. Whenever an LCP rises above `plen`,
   `vlcp[plen+1..LCP]` is written. `eorun` only lowers `plen` to the LCP of an
   entry already passed in this group. So every `vlcp[plen]` read was written
   while processing the current group, and the previous group (which differs
   between ranks and the global run) has no influence.
5. *Payload.* Seeds carry (contig, position, strand) from the payload. Routing
   by `RankForPair` and the later per-pair sort make the output independent of
   which rank emitted a seed (this is what the gold matches at 1–64 ranks
   test).

### 9.2 Lemma

Let K be the multiset of all extracted records, sorted by 40-mer, and let
`own : K → {0..P-1}` assign records to ranks. Suppose

(A) *group integrity:* any two records whose k-mers share their first 12
bases have the same owner; and
(B) *local order:* each rank holds its records sorted by k-mer, and every
record's LCP byte equals the number of leading bases its k-mer shares with
its predecessor in the rank's sequence, except that the first record of the
rank may hold any value ≤ 11.

Then the multiset of seed records produced by running FastGA's adaptive seed
merge independently on every rank, with no communication, equals the multiset
FastGA's merge produces on the whole sorted index (FREQ, masking and
SOFT_MASK fixed).

### 9.3 Proof

FastGA's merge on the whole index processes the 12-base groups of K one at a
time (facts 1–2), and its output is the union of per-group outputs. For a
group G, the output depends only on the ordered entries of G and their LCP
bytes, read as in fact 4.

Fix a group G. By (A) all of G is on one rank r. By (B) G appears contiguously
in r's sequence (the slice is sorted, and every record between two members of
G shares their 12-base prefix, so it is in G). It is in sorted order, so the
adapter on rank r forms exactly G as one prefix group (fact 1). Its entries 2..n
have the same predecessors as in the global order (their predecessors are in
G), so their LCP bytes are the global ones, all ≥ 12. Its first entry's
predecessor differs in the first 12 bases both globally and locally, so both
LCPs are ≤ 11. If G starts rank r's sequence, (B) allows any value ≤ 11, and
the code writes 0 or the repaired cross-rank value. By fact 4 the merge's
behaviour is identical for any first-entry LCP ≤ 11. Hence rank r's merge
over G reads the same entries, in the same order, with the same effective LCP
values as the global merge, and emits the same seed records, including the
frequency decisions of fact 3, which count over the entries of G and hence over
all occurrences of each LCP run. No other rank holds any record of G, so no
other rank emits seeds involving G. Summing over groups gives equality of the
seed multisets. Generating them used only rank-local data, hence no
communication. ∎

### 9.4 Applying it to the implementation

- *5-base buckets imply (A):* two k-mers sharing 12 bases share their first 5,
  so the same bucket, so the same owner (`select[bucket]`). Whole buckets per
  rank is **sufficient**; whole 12-base groups is the **exact** requirement.
- *(B):* `sort_records` + `recalc_all_lcps` give sorted order and semantic
  LCPs, with byte 0 of the first record set to 0. `run_stage2` then overwrites
  it with the LCP against the last k-mer of the nearest non-empty preceding
  rank. That value is the true global LCP and is ≤ 11 whenever the rank starts
  a new group, which (A) guarantees. The repair is not needed for seed
  equality (any value ≤ 11 is equivalent); it makes the distributed index
  byte-identical to the global one, which the M-03 tests check.
- *Contiguous ownership* is not required by the lemma. It is required for the
  global-order cross-rank check and for the concatenation of slices to *be*
  the global sorted index, which is how we test (B).
- Assumptions to state in the paper: deterministic FREQ/masking parameters
  identical on all ranks; the merge runs with group-at-a-time caching (true
  for FastGA's merge, `NTHREADS` irrelevant); masks are 0 (soft masking not
  supported, `mlen = KMER+1`).

### 9.5 Trade-off

Keeping whole buckets makes the partition coarse: with 5-base buckets the
largest bucket (AAAAA, 0.88% of the records) bounds the best possible
max/mean at 0.0088·P, which exceeds 1 from P ≈ 114. The lemma allows any
split at a 12-base group boundary, so a deeper contiguous prefix (6–8 bases)
removes the imbalance at no cost in communication. The indivisible unit then
becomes the largest 12-base group (0.175%), which would bind only beyond
~570 ranks.

### 9.6 LaTeX

```latex
\begin{lemma}[Communication-free seed generation]\label{lem:n1}
Let the extracted records be distributed over $P$ ranks such that
(i) all records whose $k$-mers share their first $12$ bases reside on the same
rank, and (ii) each rank holds its records sorted by $k$-mer, with every
record's LCP equal to the number of leading bases shared with its predecessor
on that rank, except that a rank's first record may carry any LCP $\le 11$.
Then running FastGA's adaptive seed merge independently on every rank, without
communication, yields exactly the multiset of seeds that FastGA's merge yields
on the global sorted index.
\end{lemma}

\begin{proof}
FastGA's merge partitions the sorted index into maximal groups sharing a
$12$-base prefix and processes one group at a time: it caches the group,
appends an end marker with LCP $11$, and emits seed pairs only between entries
of the cached group. Whether a seed is emitted at length $\ell \ge 12$ depends
on the number of entries of the group whose LCP run reaches $\ell$, compared
with the frequency cutoff, and on the order and LCP values of the group's
entries. The LCP of a group's first entry is used only as a lower loop bound
and is equivalent for every value $\le 11$; every other value the merge
consults is written while processing the same group. By (i), every group lies entirely
on one rank. By (ii), it appears there contiguously, in the same order, and
with the same LCPs for all but its first entry, whose LCP is $\le 11$ both
globally and locally. Hence each rank's merge processes each of its groups
exactly as the global merge does, including the complete occurrence counts for
the frequency cutoff, and no group is processed on two ranks. The seed
multisets are therefore equal, and their computation uses only rank-local
data.
\end{proof}

\noindent Partitioning by the first five bases satisfies (i), since $k$-mers
sharing $12$ bases share their first five. Local sorting followed by LCP
recomputation and the cross-rank LCP repair satisfies (ii). The price of (i)
is granularity: a bucket cannot be split, so the largest bucket (poly-A,
$0.88\%$ of the records on hsapiens90.chr6) limits balance once $P$ approaches
$1/0.0088 \approx 114$. Any finer prefix of at most $12$ bases also satisfies
(i) and removes this limit.
```

---

## 10. Answers to the R3 extraction/sort questions

1. **Do five-base buckets cause imbalance?** The buckets themselves are very
   uneven (max/mean 8.97). The per-rank slices are not, up to 64 ranks
   (1.05 / 1.18), because each rank holds many buckets. The 64-rank excess is
   the greedy splitter's (optimal split of the same buckets: 1.05). From about
   128 ranks the poly-A bucket alone exceeds a rank's share, and five-base
   partitioning cannot balance any more; a deeper prefix (≤ 12 bases) can, with
   the same correctness argument (§9).
2. **Why is sorting speedup weak from 32 to 64 ranks?** Both runs use the
   same 4 nodes. The exchange is bound by node network bandwidth and does not
   shrink (4.6 → 4.2 s). The local sort (56% of the stage) shrank only 1.31×.
   `msd_sort` sorts each part on one thread; with 256 first-byte parts a
   64-rank slice has ~5 parts, and the rank holding the poly-A bucket sorts it
   on a single thread. Giving it FastGA's own 1,024 parts and splitting heavy
   parts cut the stage by 22–29% (human 64: 17.97 → 13.13 s) with identical
   output. Efficiency was 66%; it is now 68% at fixed nodes and 79% when nodes
   double (64 r / 4 n → 128 r / 8 n).
3. **Which radix sort?** FastGA's own MSD radix sort (`MSDsort.c`): in place,
   most significant byte first, keys = the 10 packed k-mer bytes, one thread
   per pre-grouped part, LCPs recomputed afterwards (§6).
4. **Why 45–50× more memory than FastGA?** That ratio is human-only (2.4–4.3×
   on primates). FastGA builds its index out of core (4.7 GB) and its
   alignment needs 15 GB on human. We keep the 206 GB index in RAM, about three
   times over during the sort, spread over 32–64 ranks, plus 2 GiB of exchange
   staging per rank. Per rank we used 14–23 GB, about FastGA's single-process
   peak. Freeing the extraction buffer after packing and exchanging with
   `MPI_Alltoallv_c` (no staging) brought the sum to 27.6× / 29.0× (human 32 /
   64) and the per-rank peak to 14.2 / 8.7 GB. Going much lower needs a
   streaming exchange (~1× the index, §5.4) or an out-of-core index as in
   FastGA.
5. **Can communication overlap computation?** Possible at three points.
   Bounded by the exchange time: ≤ 2.2% of end-to-end on human, < 0.1% on
   primates14. Not worth an asynchronous design (§8.1).
6. **Hashing / heavy k-mers as in HySortK?** Hashing whole k-mers would split
   12-base groups and break the frequency cutoff (it would count partial
   groups), and would break the sorted index. Hashing on the 12-base prefix
   would be exact (§9.4) but gains nothing over a deeper contiguous prefix,
   which keeps the global order. Heavy handling across ranks is unnecessary:
   the largest indivisible unit (a 12-base group) is 0.175% of the records.
   Heavy buckets *within a rank* did matter. The poly-A bucket was sorted by
   one thread, and splitting it (§6) is our counterpart of HySortK's
   heavy-k-mer treatment, done without changing the order.

---

## 11. Suggested paper wording and figures

- Figure (`r3_bucket_hist.pdf`): "Records per five-base bucket,
  hsapiens90.chr6. Buckets differ by up to 85,000× (CpG-depleted vs. poly-A),
  but each rank holds many buckets: the dashed lines show a rank's ideal load
  at 64 and 128 ranks."
- Figure (`r3_balance_vs_ranks.pdf`): max/mean per-rank load vs. ranks for
  the production splitter, optimal contiguous five-base split, and eight-base
  split; extraction split alongside.
- Figure (`r3_stage2_breakdown.pdf`): sort-stage decomposition.
- Figure (`r3_sort_before_after.pdf`): sort stage and per-rank stage-2 memory,
  pre-R3 vs. R3 defaults, human 32–256 ranks and primates14 16/64 ranks.
- Text: "At 64 ranks the largest sorted slice is 18% above the mean; an
  optimal split of the same buckets reduces this to 5%. Imbalance is therefore
  not caused by the prefix distribution below about 100 ranks, but the
  largest bucket (0.88% of records) bounds five-base partitioning at higher
  rank counts; any prefix of up to 12 bases preserves Lemma~\ref{lem:n1}."
- Memory: report per-rank peak and aggregate side by side, and FastGA's
  index and alignment peaks separately. State that the aggregate counts the
  in-memory index replicas and that FastGA's index is out of core.
- Sort stage, from §7.3: "Freeing the record buffer after packing,
  exchanging with a single large-count \texttt{MPI\_Alltoallv\_c}, giving the
  radix sort FastGA's 1,024-part layout and splitting heavy (poly-A/T) parts
  reduce the sort stage by 22--29\% and its per-rank memory by 38\% (hsapiens90,
  64 ranks: 17.97 to 13.13 s, 14.1 to 8.7 GB), with byte-identical output."
- End-to-end: report the sort-stage and memory gains; the end-to-end time is
  alignment-bound, and these gains are within the alignment's run-to-run
  spread.

---

## 12. P1-C extraction balance and P1-D scaling beyond 127 ranks

### 12.1 Extraction (whole-contig split)

| | 32 ranks | 64 ranks | 128 | 256 |
|---|---|---|---|---|
| human: records max/mean (idle ranks) | 1.12 (0) | 1.35 (0) | 1.75 (0) | 2.84 (29) |
| primates14 | 2.96 (10) | 4.59 (36) | 9.17 (98) | 18.3 (223) |
| primates4 | 4.28 (11) | 8.55 (33) | 17.1 (85) | 34.2 (206) |

The cause is contig granularity (human: 1,421 contigs; primates14: 81 with
the largest 4.6× a 64-rank share). Fallback reads are negligible: 260 of
195.7 M records on human rank 0 needed a `Get_Contig_Piece` outside the scan
window. Cost today: the wait for the slowest extractor is 1.0 / 1.5 s
(human 32 / 64) and 2.3–2.7 s (primates14), ≤ 0.8% of end-to-end on human
and ≤ 0.2% on primates14. It also feeds pack time on the heaviest extractor.
**No change** at ≤ 64 ranks. For ≥ 128 ranks, splitting contigs at base
offsets is required. That is an extraction rewrite, because a split must
re-create the syncmer window state (4 × 12-mers) and the k-mer carry at the
cut; the existing multi-window path (`SCAN_MAX` windows with `KMER_CARRY`)
already does this within one rank and is the template.

### 12.2 What limits extraction/sort beyond 127 ranks

- **Alignment side (not touched):** `RankForPair` is `int8_t` (≤ 127 ranks)
  and two `ncontig²` maps are live (`contig_assignment.c`,
  `pga-mpi-merge.c`). Saksham's area.
- **Sort side, no hard cap:** `compute_destination_map`, the Allgather in
  the boundary check, and the exchange all work for any P; the bench ran 128
  and 256 ranks (4 and 8 nodes) with identical digests.
- **Soft limits, in order of onset (human):** whole-contig extraction (1.35 at
  64, 1.75 at 128); five-base buckets (floor 1.12 at 128, 2.24 at 256; P >
  1,024 leaves ranks empty). Measured: 256 ranks on 8 nodes is slower than
  128 (8.79 vs 8.34 s). The slowest rank still holds AAAAA, whose largest
  9-base group (AAAAAAAAA, 33 M records) is now its indivisible sort unit.
  Per-rank fixed memory remains (128 MiB adapter index, MPI).
- **Needed for 128+ ranks on the extraction/sort side:**
  1. deeper contiguous prefix buckets (8 bases = 65,536 buckets: the
     histogram Allreduce grows to 512 KiB, everything else is unchanged); the
     sort keeps its own ≤ 1,024 parts;
  2. a sub-contig extraction split;
  3. sort threads = physical cores per rank, set in the job script.

  (1) and (2) are not implemented. Their benefit is nil below ~100 ranks, and
  their correctness rests on §9 (for 1) and the existing window logic (for 2).
  The large-count exchange and the heavy-part split are done.

---

## 13. Risks and next steps

1. **Merge:** open a PR from `r3-extract-sort`. It touches `pga-mpi-merge.c`
   only in the extraction/sort block (~lines 330–375); tell Saksham before he
   rebases O3.
2. **Report numbers:** the report's stage tables come from the pre-R3 sweep.
   Re-running the 12-configuration sweep on this branch would give
   end-to-end numbers with the new defaults (sort and memory change;
   alignment does not).
3. **Job scripts:** set `PGA_SORT_THREADS` to the rank's physical cores
   (`-c / 2`); measured gain 0–12% of the sort stage on human/primates14 (18% on yeast7), gold-validated (`phys`).
4. **Deferred** (only for ≥ 128 ranks): deeper prefix buckets (§12.2), a
   sub-contig extraction split, a second level of heavy-part splitting for
   AAAAAAAAA. Streaming exchange only if memory becomes binding.
5. **Upstream:** report the `msd_sort` threading defects (§6) to FastGA's
   authors.
6. **Risks:** the order of equal k-mers changes with the new sort layout.
   That order was never specified, and 32/32 gold runs and the O1-0 shuffle
   test show the output does not depend on it. Any future code that reads
   equal-k-mer order would need re-validation. `MPI_Alltoallv_c` needs
   MPI ≥ 4 (compile-time guard; older MPIs fall back to rounds).
