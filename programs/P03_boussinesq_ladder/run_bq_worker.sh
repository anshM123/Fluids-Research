#!/bin/bash
# 2D between-rung spectra, remaining states split over two workers (states re-converged at their exact λ)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for lam in "$@"; do
  f=$(grep "^C " exact_lam.txt | awk -v l=$lam '$3==l {print $2}')
  python3 reconverge.py bq $f $lam 0.025 Y2_Cc >> reconverge.out 2>&1
  python3 bq_flow_spec.py Y2_Cc_lam$lam.npy $lam 0.025 3.6 0.1 -30 0.2 >> bq_flow.out 2>&1
done
