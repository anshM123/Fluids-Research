#!/bin/bash
# Hou–Luo between-rung spectra on states re-converged at their exact λ (E states, z = 10.3–13.1)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
grep "^E " exact_lam.txt | while read k f lam; do
  python3 reconverge.py hl $f $lam 32768 40 160 -30 hl_q_Ec >> reconverge.out 2>&1
  HLGRID="32768,40,160" python3 hl_deep_spec.py -60 3.2 0.05 hl_q_Ec_lam$lam.npy >> hl_saw_exact.out 2>&1
done
