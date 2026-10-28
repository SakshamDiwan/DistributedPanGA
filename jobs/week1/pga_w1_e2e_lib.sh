# Shared by Job C and D. Runs pga-mpi-merge at N ranks in a fresh dir and checks
# the full sorted (duplicate-preserving) aggregate PAF against FastGA gold.
# Usage: run_one <N> <dir> ; needs $R $STEM $GOLD_SORTED set.
run_one() {
  local N=$1 D=$2
  rm -rf "$D"; mkdir -p "$D"; cd "$D" || return 1
  local t0=$(date +%s.%N)
  srun -n "$N" -c 16 --cpu-bind=verbose,cores "$HOME/pga-jobs/week1/pga_w1_rankwrap.sh" "$R/pga-mpi-merge" "$STEM" > merge.out 2> merge.err
  local rc=$?
  local t1=$(date +%s.%N)
  local wall=$(echo "$t1 - $t0" | bc)
  cat out.rank*.paf 2>/dev/null | LC_ALL=C sort > got.sorted.paf
  local n=$(wc -l < got.sorted.paf) dup=$(uniq -d got.sorted.paf | wc -l)
  local empty=0; for f in out.rank*.paf; do [ -s "$f" ] || empty=$((empty+1)); done
  local nfiles=$(ls out.rank*.paf 2>/dev/null | wc -l)
  local conv=$(grep -c 'ALNtoPAF failed' merge.err)
  local verdict
  if [ "$n" -gt 0 ] && cmp -s got.sorted.paf "$GOLD_SORTED"; then verdict=MATCH; else verdict=DIFFER; fi
  local maxrss=$(grep -h 'Maximum resident' time.rank* 2>/dev/null | awk '{print $NF}' | sort -n | tail -1)
  echo "RESULT N=$N rc=$rc wall_s=$wall alignments=$n dup_lines=$dup rank_files=$nfiles empty_rank_files=$empty alntopaf_fail=$conv max_rank_rss_kb=$maxrss gold=$verdict dir=$D"
  echo "CSV    $(tail -1 pga-mpi-merge.timing.csv 2>/dev/null)"
  grep -hE 'own slice has|adapter .* entries|received N=' merge.err | sed 's/^/  work: /'
  for f in time.rank*; do echo "  rss: $f $(grep 'Maximum resident' $f | awk '{print $NF}') kB  elapsed $(grep 'Elapsed (wall' $f | awk '{print $NF}')"; done
  if [ "$verdict" = DIFFER ]; then diff got.sorted.paf "$GOLD_SORTED" | head -8 | sed 's/^/  diff: /'; fi
  cd - >/dev/null
}
