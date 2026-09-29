#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 glob_run.py Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.252349.npy 1.2523487 l2_h025_Nb24 0.025 24 0.025 -8 8,10,12 > /dev/null 2>&1
python3 glob_run.py Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.184253.npy 1.1842533 l3_h025_Nb24 0.025 24 0.025 -8 8,10,12 > /dev/null 2>&1
python3 glob_run.py Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.144986.npy 1.1449857 l4_h025_Nb24 0.025 24 0.0125 -8 8,10,12 > /dev/null 2>&1
