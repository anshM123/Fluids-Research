#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 stab_contour2.py Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.101582.npy 1.1015816513 0.0125 C6 0.065 1.2 0.8 12 -20 > /dev/null 2>&1
