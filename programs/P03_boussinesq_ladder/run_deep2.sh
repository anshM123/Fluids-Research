#!/bin/bash
# deep branch states with 28 eigenvalues of T_μ per point (the top-16 set saturates for z ≳ 15)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 HLK=28 HLNUMIN=0.5
log=$1; shift
for lam in "$@"; do
  python3 hl_deep_spec.py -100 4.5 0.1 hl_q_F_lam$lam.npy >> $log 2>&1
done
