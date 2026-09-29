#!/bin/bash
# grid variant of the λ₇ right-half-plane count (hs = 0.0125), after the low-μ Hou–Luo probe
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while kill -0 4031 2>/dev/null; do sleep 30; done
python3 stab_contour2.py Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.088338.npy 1.0883384 0.0125 C7 0.055 1.2 0.8 12 -30 > /dev/null 2>&1
