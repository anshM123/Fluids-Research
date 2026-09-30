#!/bin/bash
# spectra at the branch points between U_2 and U_3 (stage 2b P5), launched as the states appear
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
for spec in "0.2927270 b025 3.6" "0.2734243 b050 3.8" "0.2565098 b075 4.1"; do
  set -- $spec
  until [ -f ipm_bp_lam$1.npy ]; do sleep 20; done
  sleep 5
  timeout 7000 python3 ipm_flow_spec.py ipm_bp_lam$1.npy $1 $3 0.15 -30 0.05 12 > ipm_spec_23$2.out 2>&1
done
