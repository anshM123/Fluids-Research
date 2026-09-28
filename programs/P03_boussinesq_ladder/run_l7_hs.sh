#!/bin/sh
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 30000 python3 bq_crossing.py Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.089008.npy 1.0890078 Y2_B_cross_lam1.0800.npy 1.08 32 0.0125 -20 100 32 0.025 > cross_l7_hs0125.log 2>&1
