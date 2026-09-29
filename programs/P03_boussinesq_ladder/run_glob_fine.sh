#!/bin/bash
# hs = 0.0125 global runs for the deep profiles (after run_glob_all.sh has finished)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
until grep -q "geometric" glob_g7_h025_Nb24_s12.log 2>/dev/null; do sleep 30; done
P="Ycross_Nb32"
python3 glob_run.py ${P}_hs0.0125_ss-20.0_sm100.0_lam1.144986.npy 1.1449857 g4_h0125_Nb24_s12 0.0125 24 0.0125 -12 8,10,12 > /dev/null 2>&1
python3 glob_run.py ${P}_hs0.0125_ss-20.0_sm100.0_lam1.119474.npy 1.1194738 g5_h0125_Nb24_s12 0.0125 24 0.0125 -12 8,10,12 > /dev/null 2>&1
python3 glob_run.py ${P}_hs0.0125_ss-20.0_sm100.0_lam1.101582.npy 1.1015817 g6_h0125_Nb24_s12 0.0125 24 0.0125 -12 8,10,12 > /dev/null 2>&1
python3 glob_run.py ${P}_hs0.0125_ss-20.0_sm100.0_lam1.088338.npy 1.0883384 g7_h0125_Nb24_s12 0.0125 24 0.0125 -12 8,10,12 > /dev/null 2>&1
python3 glob_run.py ${P}_hs0.025_ss-20.0_sm100.0_lam1.184253.npy 1.1842533 g3_h0125_Nb24_s12 0.0125 24 0.025 -12 8,10,12 > /dev/null 2>&1
