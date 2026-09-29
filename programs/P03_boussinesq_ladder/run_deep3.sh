#!/bin/bash
while kill -0 4013 2>/dev/null; do sleep 30; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 HLK=28 HLNUMIN=0.3
python3 hl_deep_spec.py -100 4.5 0.1 hl_q_F_lam1.064267.npy >> hl_deep_E.out 2>&1
