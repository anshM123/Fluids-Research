#!/bin/bash
# 2D spectral flow between rungs 3–4 and 4–5 (hs = 0.025 branch states of scan C), after the C6 contour finishes.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while kill -0 2338 2>/dev/null; do sleep 30; done
for s in "Y2_C_lam1.1576.npy 1.1576" "Y2_C_lam1.1684.npy 1.1684" "Y2_C_lam1.1322.npy 1.1322" "Y2_C_lam1.1481.npy 1.1481" "Y2_C_lam1.1397.npy 1.1397"; do
  set -- $s
  python3 bq_flow_spec.py $1 $2 0.025 3.6 0.1 -30 0.2 >> bq_flow.out 2>&1
done
