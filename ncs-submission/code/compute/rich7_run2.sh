#!/bin/bash
# Richardson / grid-shift queue for λ7: 3 concurrent solves; jobs needing an h7 state wait until it exists.
cd "$(dirname "$0")"; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run_one() {
  line="$1"; w=$(echo "$line" | sed -n 's/^WAIT=\([^ ]*\) .*/\1/p')
  [ -n "$w" ] && until [ -f "$w" ]; do sleep 30; done
  cmd=$(echo "$line" | sed 's/^WAIT=[^ ]* //')
  tag=$(echo "$cmd" | awk '{print $(NF-2)"_d"$(NF-1)"_hs"$NF}')
  [ -s "rich7_$tag.out" ] && grep -q "m-2=" "rich7_$tag.out" && exit 0     # done before a restart
  eval "$cmd" > "rich7_$tag.out" 2>&1
}
export -f run_one
grep -v '^#' rich7_queue2.txt | xargs -d '\n' -P 3 -I{} bash -c 'run_one "$@"' _ {}
