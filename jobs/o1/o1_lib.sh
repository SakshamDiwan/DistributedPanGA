# Like week1 run_one, but with an explicit binary and GDB stem.
# Usage: run_bin <label> <binary> <stem> <N> <dir> ; needs $GOLD_SORTED.
run_bin() {
  local label=$1 BIN=$2 STM=$3 N=$4 D=$5
  rm -rf "$D"; mkdir -p "$D"; cd "$D" || return 1
  local t0=$(date +%s.%N)
  srun -t ${RUN_LIMIT_MIN:-10} -n "$N" -c 16 --cpu-bind=cores "$HOME/pga-jobs/week1/pga_w1_rankwrap.sh" "$BIN" "$STM" > merge.out 2> merge.err
  local rc=$?
  local wall=$(echo "$(date +%s.%N) - $t0" | bc)
  cat out.rank*.paf 2>/dev/null | LC_ALL=C sort > got.sorted.paf
  local n=$(wc -l < got.sorted.paf) dup=$(uniq -d got.sorted.paf | wc -l)
  local conv=$(grep -c 'ALNtoPAF failed' merge.err)
  local verdict=DIFFER; [ "$n" -gt 0 ] && cmp -s got.sorted.paf "$GOLD_SORTED" && verdict=MATCH
  local sys=$(awk '/System time/{s+=$NF} END{printf "%.1f", s}' time.rank*)
  local usr=$(awk '/User time/{s+=$NF} END{printf "%.1f", s}' time.rank*)
  local maxrss=$(grep -h 'Maximum resident' time.rank* | awk '{print $NF}' | sort -n | tail -1)
  echo "RESULT $label N=$N rc=$rc wall_s=$wall alignments=$n dup_lines=$dup alntopaf_fail=$conv user_s_all=$usr sys_s_all=$sys max_rank_rss_kb=$maxrss gold=$verdict dir=$D"
  echo "CSV    $label $(tail -1 pga-mpi-merge.timing.csv 2>/dev/null)"
  o1_summary "$D"
  [ "$verdict" = DIFFER ] && diff got.sorted.paf "$GOLD_SORTED" | head -6 | sed 's/^/  diff: /'
  cd - >/dev/null
}
