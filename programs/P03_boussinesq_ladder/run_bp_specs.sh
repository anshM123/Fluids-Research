#!/bin/bash
# launch the between-rung spectra as soon as the states exist (P5 test of PREDICTIONS_IPM.md stage 2a)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
for spec in "0.3778542 050 2.8" "0.3435536 075 3.1"; do
  set -- $spec
  until [ -f ipm_bp_lam$1.npy ]; do sleep 15; done
  sleep 5
  timeout 7000 python3 ipm_flow_spec.py ipm_bp_lam$1.npy $1 $3 0.1 -30 0.05 12 > ipm_spec_bp$2.out 2>&1 &
done
wait
