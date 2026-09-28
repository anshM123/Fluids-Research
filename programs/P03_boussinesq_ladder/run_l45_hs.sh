#!/bin/sh
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 30000 python3 bq_crossing.py Y2_B_lam1.1200.npy 1.12 Y2_B_cross_lam1.1100.npy 1.11 32 0.0125 -20 100 32 0.025 > cross_l5_hs0125.log 2>&1 &
timeout 30000 python3 bq_crossing.py Y2_B_lam1.1600.npy 1.16 Y2_B_cross_lam1.1400.npy 1.14 32 0.0125 -20 100 32 0.025 > cross_l4_hs0125.log 2>&1 &
wait
