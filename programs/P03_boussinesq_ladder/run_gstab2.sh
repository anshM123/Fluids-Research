#!/bin/bash
# λ7 stability on its converged hs = 0.0125 state, then the s_min = −20 repeats for all profiles
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
until grep -q RESULT gstab_g6.log 2>/dev/null; do sleep 20; done
SH="0.1:0,0.2:0,0.35:0,0.5:0,0.65:0,0.8:0,1.0:0,0.2:0.4,0.5:0.5,0.9:0.5,0.3:1.0,0.8:1.0"
SH7="0.07:0,0.1:0,0.2:0,0.35:0,0.5:0,0.65:0,0.8:0,1.0:0,0.2:0.4,0.5:0.5,0.9:0.5,0.3:1.0,0.8:1.0"
python3 glob_stab.py glob_g7_h0125_Nb24_s12_smax12.npy -12 12 0.0125 24 g7 $SH7 12 > /dev/null 2>&1
F=(Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.920559.npy Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.399096.npy Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.252349.npy Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.184253.npy Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.144986.npy Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.119474.npy Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.101582.npy Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam1.088338.npy)
L=(1.9205593 1.3990961 1.2523487 1.1842533 1.1449857 1.1194738 1.1015817 1.0883384)
SRC=(0.025 0.025 0.025 0.025 0.0125 0.0125 0.0125 0.0125)
HG=(0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.0125)
for n in 1 2 3 4 5 6 7 0; do
  python3 glob_run.py ${F[$n]} ${L[$n]} g${n}_s20 ${HG[$n]} 24 ${SRC[$n]} -20 12 > /dev/null 2>&1
  S=$SH; [ $n = 7 ] && S=$SH7
  python3 glob_stab.py glob_g${n}_s20_smax12.npy -20 12 ${HG[$n]} 24 g${n}s20 $S 12 > /dev/null 2>&1
done
