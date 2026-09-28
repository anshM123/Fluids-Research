#!/bin/sh
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
until grep -q CROSSING cross_l5_hs0125.log 2>/dev/null && grep -q CROSSING cross_l4_hs0125.log 2>/dev/null; do sleep 30; done
timeout 40000 python3 bq_crossing.py Y2_B_cross_lam1.1100.npy 1.1025 Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.101523.npy 1.1010 32 0.00625 -20 100 32 0.025 > cross_l6_hs00625.log 2>&1
