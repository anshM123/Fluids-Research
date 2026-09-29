#!/bin/bash
# independent reproduction of λ0–λ7 with the global solver: s_min = −12, s_max = 8, 10, 12 (+ extrapolation)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
P="Ycross_Nb32"
python3 glob_run.py ${P}_hs0.025_ss-20.0_sm100.0_lam1.920559.npy 1.9205593 g0_h025_Nb24_s12 0.025 24 0.025 -12 8,10,12 > /dev/null 2>&1
python3 glob_run.py ${P}_hs0.025_ss-20.0_sm100.0_lam1.399096.npy 1.3990961 g1_h025_Nb24_s12 0.025 24 0.025 -12 8,10,12 > /dev/null 2>&1
python3 glob_run.py ${P}_hs0.025_ss-20.0_sm100.0_lam1.252349.npy 1.2523487 g2_h025_Nb24_s12 0.025 24 0.025 -12 8,10,12 > /dev/null 2>&1
python3 glob_run.py ${P}_hs0.025_ss-20.0_sm100.0_lam1.184253.npy 1.1842533 g3_h025_Nb24_s12 0.025 24 0.025 -12 8,10,12 > /dev/null 2>&1
python3 glob_run.py ${P}_hs0.0125_ss-20.0_sm100.0_lam1.144986.npy 1.1449857 g4_h025_Nb24_s12 0.025 24 0.0125 -12 8,10,12 > /dev/null 2>&1
python3 glob_run.py ${P}_hs0.0125_ss-20.0_sm100.0_lam1.119474.npy 1.1194738 g5_h025_Nb24_s12 0.025 24 0.0125 -12 8,10,12 > /dev/null 2>&1
python3 glob_run.py ${P}_hs0.0125_ss-20.0_sm100.0_lam1.101582.npy 1.1015817 g6_h025_Nb24_s12 0.025 24 0.0125 -12 8,10,12 > /dev/null 2>&1
python3 glob_run.py ${P}_hs0.0125_ss-20.0_sm100.0_lam1.088338.npy 1.0883384 g7_h025_Nb24_s12 0.025 24 0.0125 -12 8,10,12 > /dev/null 2>&1
