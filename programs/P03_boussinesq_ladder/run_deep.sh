#!/bin/bash
# usage: run_deep.sh WAIT_PID LOG STATES...   (waits for WAIT_PID to exit, then scans the states one by one)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
w=$1; log=$2; shift 2
while kill -0 $w 2>/dev/null; do sleep 30; done
for lam in "$@"; do
  python3 hl_deep_spec.py -100 4.4 0.1 hl_q_F_lam$lam.npy >> $log 2>&1
done
