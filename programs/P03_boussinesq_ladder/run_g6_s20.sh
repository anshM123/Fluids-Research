#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
until grep -q RESULT gstab_g7b_s20.log 2>/dev/null; do sleep 20; done
python3 glob_run.py Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.101582.npy 1.1015817 g6b_s20 0.0125 24 0.0125 -20 10,12 > /dev/null 2>&1
python3 glob_run.py Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.119474.npy 1.1194738 g5b_s20 0.0125 24 0.0125 -20 12 > /dev/null 2>&1
