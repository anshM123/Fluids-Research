#!/bin/bash
export OMP_NUM_THREADS=1
# domain variant (η0 = −40) and contour variant (x_lo = 0.04) for rungs 5–10
for f in hl_cross_E32768_lam1.13772677.npy hl_cross_E32768_lam1.11731019.npy hl_cross_E32768_lam1.10214553.npy hl_cross_E32768_lam1.09044036.npy hl_cross_E32768_lam1.08113374.npy hl_cross_E32768_lam1.07355524.npy; do
  python3 hl_contour.py $f 32768 0.04 1.5 1.5 20 -40 2>&1 | grep "HL λ"
done
# grid variant N = 65536 for rungs 8 and 10
for f in hl_cross_E32768_lam1.09044036.npy hl_cross_E32768_lam1.07355524.npy; do
  python3 hl_refine_N.py $f 32768 65536 2>&1 | grep "λ="
done
for f in hl_crossN65536_lam1.09044036.npy hl_crossN65536_lam1.07355524.npy; do
  python3 hl_contour.py $f 65536 0.03 1.5 1.5 20 -30 2>&1 | grep "HL λ"
done
