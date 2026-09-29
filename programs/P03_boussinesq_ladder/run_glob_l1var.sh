#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
F=Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.399096.npy
python3 glob_run.py $F 1.3990961 l1_h0125_Nb24 0.0125 24 0.025 -8 8,10,12 > /dev/null 2>&1
python3 glob_run.py $F 1.3990961 l1_h025_Nb32 0.025 32 0.025 -8 8,10,12 > /dev/null 2>&1
python3 glob_run.py $F 1.3990961 l1_h025_Nb24_smin10 0.025 24 0.025 -10 8,10,12 > /dev/null 2>&1
