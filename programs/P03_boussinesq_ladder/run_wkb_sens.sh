#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 wkb_sens.py Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.144986.npy Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.144986.npy Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.119474.npy Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.101582.npy Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.088338.npy Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.184253.npy > wkb_sens.out 2>&1
