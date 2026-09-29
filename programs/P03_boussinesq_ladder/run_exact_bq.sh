#!/bin/bash
# 2D between-rung spectra on states re-converged at their exact λ (scan C, hs = 0.025)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
grep "^C " exact_lam.txt | while read k f lam; do
  python3 reconverge.py bq $f $lam 0.025 Y2_Cc >> reconverge.out 2>&1
  python3 bq_flow_spec.py Y2_Cc_lam$lam.npy $lam 0.025 3.6 0.1 -30 0.2 >> bq_flow.out 2>&1
done
