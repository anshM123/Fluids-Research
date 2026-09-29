#!/bin/bash
# Hou–Luo spectrum at branch points between the rungs n = 7…10 (E states, N = 32768), to see how the ladder offset
# behaves away from the smooth profiles.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 HLGRID="32768,40,160"
for lam in 1.096993 1.093371 1.090009 1.086881 1.083963 1.081235 1.078678 1.076278; do
  nice -n 5 python3 hl_deep_spec.py -60 3.2 0.05 hl_q_E_lam$lam.npy >> hl_saw.out 2>&1
done
