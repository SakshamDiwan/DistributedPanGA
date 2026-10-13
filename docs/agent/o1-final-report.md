# O1 final report — HySortK-inspired extraction (closed)

An ACM two-column version is in `docs/o1-report/o1-final-report.tex` (built PDF:
`docs/o1-report/o1-final-report.pdf`; rebuild with
`module load texlive/2024; latexmk -pdf o1-final-report.tex`).

**Status: closed, 2026-10-12.** O1 was the proposal's first objective: speed
up the extraction stage of DistributedPanGA, using ideas from HySortK (Li &
Guidi, ICPP 2024). The deliverable is merged in PR #7. The other planned
sub-steps were dropped because measurements showed they could not pay off.
Work that started from O1 but belongs elsewhere is listed at the end.

## 1. What was done

| step | what | result (yeast7, 1 node) | where |
|---|---|---|---|
| Week 1 | Verified the pipeline against independent FastGA gold (43,126 alignments); baseline | all phase checks and end-to-end runs match gold | PR #7, `baseline.md` |
| O1-0 | Per-rank instrumentation and a pre-sort shuffle test | found one `lseek` system call per extracted record; work was balanced; output does not depend on pre-sort order | PR #7, `o1-0-results.md` |
| **O1a** | Pack each seed's k-mer from the sequence already in memory instead of re-reading it from the file | **extraction 49.8 s → 1.5 s (N=1), 28.6 s → 0.2 s (N=8): 34–136× faster; whole run 1.7–3.7× faster**; 1→8-rank speedup 2.9× → 6.2× | PR #7, `o1a-results.md` |

Correctness: every O1 run matched FastGA gold (20 of 20 O1a runs). The
extraction output is byte-identical to FastGA's own index.

## 2. What was not done, and why

| planned step | decision | evidence |
|---|---|---|
| O1b — multithreaded extraction | dropped | After O1a, extraction is ~2% of run time (0.2 s of 10.6 s at N=8). Threads could save at most ~1.5 s at N=1. |
| O1c — destination-aware staging (HySortK's per-thread × per-destination buffers) | dropped | The speed gain is ≈0.45 s. When proposed, it could not lower peak memory either: the peak was in the alignment stage. After PR #9 freed those buffers, stage 2 is the peak at 1 rank (2.96 GiB), so O1c could now save some memory. That is optional and outside O1's speed goal. |
| O1d — batched, overlapped exchange (HySortK's `SendingGroup`/`ReceivingGroup`) | dropped | The k-mer exchange takes ≤ 0.13 s on one node, so there is nothing to hide. Very large exchanges were solved separately by bounded rounds (PR #10). |
| Supermer compression | not applicable | Seeds are routed by their own first bases, and consecutive selected seeds go to the same rank with probability ≈ 1/P. Supermers would not compress, and routing by minimizer would break correctness. (`o1-design.md` §6) |

## 3. Conclusion on HySortK

- **Did not transfer:**
  - Supermers and heavy-hitter handling (k-mer counting assumptions).
  - Extraction-side staging and overlap: the measured extraction cost was a
    per-record file seek, not something HySortK addresses.
  - HySortK's **static, size-based task dispatch**, tested later on the
    alignment stage. Alignment cost follows how much sequence is aligned,
    not the size of the input, so it was no better than the existing
    assignment (`decision-gate-1.md`, `decision-gate-2.md`).
- **Did transfer:** HySortK's **task abstraction**, i.e. many independent
  tasks, more than ranks. It frames the follow-on scheduling work.

## 4. Follow-on work (not O1)

| item | status |
|---|---|
| `ALNtoPAF` timeout / failure handling | merged, PR #8 |
| Merge/alignment instrumentation and early memory frees (−25 to −38% peak) | PR #9 (open) |
| MPI exchanges larger than 2 GB per rank | PR #10 (open); needed for every primate dataset |
| Scratch files off `/tmp`; per-task cost diagnostics; decision gates | PR #11 (open) |
| Larger datasets with FastGA gold: scerevisiae8, primates4.chr6, primates14.chr6 (incl. a 2-node run) | all match gold |
| **Alignment-stage scheduling (proposal O3)** | Evidence gathered; not started. A dynamic queue with self-alignments of long contigs first models at −33 to −42% for the slowest rank at 16–32 ranks. |
| Splitting large self-alignment tasks | Open research question; needed to scale primates14 past ~16 ranks. |
| 127-rank limit (8-bit contig-pair map) | Not fixed; needed above 8 nodes × 16 ranks. |
