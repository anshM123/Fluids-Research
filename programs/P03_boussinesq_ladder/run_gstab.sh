#!/bin/bash
# independent stability (global discretisation) for λ0–λ7; real and complex shifts
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
SH="0.1:0,0.2:0,0.35:0,0.5:0,0.65:0,0.8:0,1.0:0,0.2:0.4,0.5:0.5,0.9:0.5,0.3:1.0,0.8:1.0"
for n in 0 1 2 3 4 5 6 7; do
  until [ -f glob_g${n}_h025_Nb24_s12_smax12.npy ]; do sleep 30; done
  python3 glob_stab.py glob_g${n}_h025_Nb24_s12_smax12.npy -12 12 0.025 24 g${n} $SH 12 > /dev/null 2>&1
done
