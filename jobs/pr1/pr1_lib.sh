# run_pr1 <label> <binary> <stem> <gold.paf> <N> <dir>
# One pga-mpi-merge run: per-rank /usr/bin/time, gold check (sorted, duplicates kept),
# and (for instrumented binaries) the alignment-balance report.
run_pr1() {
  local label=$1 BIN=$2 STM=$3 GOLD=$4 N=$5 D=$6
  rm -rf "$D"; mkdir -p "$D"; cd "$D" || return 1
  local t0=$(date +%s.%N)
  srun -t ${RUN_LIMIT_MIN:-15} -n "$N" -c 16 --cpu-bind=cores "$HOME/pga-jobs/week1/pga_w1_rankwrap.sh" "$BIN" "$STM" > merge.out 2> merge.err
  local rc=$?
  local wall=$(echo "$(date +%s.%N) - $t0" | bc)
  cat out.rank*.paf 2>/dev/null | LC_ALL=C sort > got.sorted.paf
  local n=$(wc -l < got.sorted.paf)
  local verdict=DIFFER; [ "$n" -gt 0 ] && cmp -s got.sorted.paf "$GOLD" && verdict=MATCH
  local maxrss=$(grep -h 'Maximum resident' time.rank* | awk '{print $NF}' | sort -n | tail -1)
  echo "RESULT $label N=$N rc=$rc wall_s=$wall alignments=$n max_rank_rss_kb=$maxrss gold=$verdict dir=$D"
  echo "CSV    $label $(tail -1 pga-mpi-merge.timing.csv 2>/dev/null)"
  if grep -q '^\[pr1\]' merge.err; then
    python3 -I $R/scripts/analyze_align_balance.py "$D" > balance.txt 2>&1; echo "  analyze rc=$?"
    grep -E '^max/mean|^per task \(|^per rank|LPT|actual max' balance.txt | sed 's/^/  /'
    grep -h '^\[pr1-mem\]' merge.err | grep -oE 'overall_peak_kb=[0-9]+' | sort -t= -k2 -n | tail -1 | sed 's/^/  max rank /'
  fi
  rm -f got.sorted.paf
  cd - >/dev/null
}
