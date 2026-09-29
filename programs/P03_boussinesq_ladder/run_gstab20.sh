#!/bin/bash
# same as run_gstab.sh with the origin truncation moved from s_min = −12 to −20 (spurious modes move, genuine ones do not)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
until grep -q RESULT gstab_g7.log 2>/dev/null; do sleep 30; done
SH="0.1:0,0.2:0,0.35:0,0.5:0,0.65:0,0.8:0,1.0:0,0.2:0.4,0.5:0.5,0.9:0.5,0.3:1.0,0.8:1.0"
F=(Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.920559.npy Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.399096.npy Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.252349.npy Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.184253.npy Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.144986.npy Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.119474.npy Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.101582.npy Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.088338.npy)
L=(1.9205593 1.3990961 1.2523487 1.1842533 1.1449857 1.1194738 1.1015817 1.0883384)
H=(0.025 0.025 0.025 0.025 0.0125 0.0125 0.0125 0.0125)
for n in 0 1 2 3 4 5 6 7; do
  python3 glob_run.py ${F[$n]} ${L[$n]} g${n}_h025_Nb24_s20 0.025 24 ${H[$n]} -20 12 > /dev/null 2>&1
  python3 glob_stab.py glob_g${n}_h025_Nb24_s20_smax12.npy -20 12 0.025 24 g${n}s20 $SH 12 > /dev/null 2>&1
done
