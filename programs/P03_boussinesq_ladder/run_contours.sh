#!/bin/bash
# argument-principle counts of all eigenvalues (real and complex) in [0.06,1.5]x[-1.5,1.5]
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 stab_contour.py Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.252349.npy 1.2523487 0.025 l2 0.06 1.5 1.5 > /dev/null 2>&1
python3 stab_contour.py Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.184253.npy 1.1842533 0.025 l3 0.06 1.5 1.5 > /dev/null 2>&1
python3 stab_contour.py Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.144986.npy 1.1449864 0.025 l4 0.06 1.5 1.5 > /dev/null 2>&1
python3 stab_contour.py Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.119482.npy 1.1194818 0.025 l5 0.06 1.5 1.5 > /dev/null 2>&1
