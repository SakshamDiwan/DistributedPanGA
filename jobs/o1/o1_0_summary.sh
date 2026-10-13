# Summarise the per-rank "[o1-0]" lines of one pga-mpi-merge run directory.
# Usage: o1_summary <dir>   (reads <dir>/merge.err)
o1_summary() {
  grep -h '^\[o1-0\]' "$1/merge.err" | sed -E 's/^\[o1-0\] rank ([0-9]+):/rank=\1/' | awk '
  { for (i = 1; i <= NF; i++) { split($i, kv, "="); v[kv[1]] = kv[2] }
    n++
    for (k in v) if (k != "rank") { s[k] += v[k]; if (n == 1 || v[k] > mx[k]) mx[k] = v[k]; if (n == 1 || v[k] < mn[k]) mn[k] = v[k] } }
  END {
    split("extract_s records pack_reads window_packs window_reads bases contigs hist_s pack_s exchange_s shuffle_s local_sort_s verify_s sort_stage_s slice_records", ks, " ")
    printf "  O1-0 ranks=%d", n
    for (i = 1; i in ks; i++) { k = ks[i]; printf " %s[min/mean/max]=%.6g/%.6g/%.6g", k, mn[k], s[k]/n, mx[k] }
    printf " pack_reads_per_record=%.4f extract_us_per_record=%.4f\n", s["pack_reads"]/s["records"], 1e6*s["extract_s"]/s["records"]
  }'
}
