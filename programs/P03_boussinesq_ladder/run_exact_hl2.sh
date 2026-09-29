#!/bin/bash
# deep Hou–Luo branch states re-converged at their exact λ (F states, z = 10.6–38)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 HLK=28
grep "^F " exact_lam.txt | while read k f lam; do
  python3 reconverge.py hl $f $lam 65536 25 75 -20 hl_q_Fc >> reconverge.out 2>&1
  HLNUMIN=0.3 python3 hl_deep_spec.py -100 4.5 0.1 hl_q_Fc_lam$lam.npy >> hl_deep_exact.out 2>&1
done
