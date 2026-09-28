#!/bin/bash
cd /home/user/Fluids-Research/programs/P02_singularity_ladders
export OMP_NUM_THREADS=1
timeout 40000 python3 mapped_continuation.py mapped_l2_hs0.02_eps0.025.npy 0.4713242277712 3000 dn -1 > mc_dn.log 2>&1 &
timeout 40000 python3 mapped_continuation.py mapped_l2_hs0.02_eps0.025.npy 0.4713242277712 400 up 1 > mc_up.log 2>&1 &
