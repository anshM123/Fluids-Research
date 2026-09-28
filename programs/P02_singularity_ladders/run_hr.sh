#!/bin/bash
cd /home/user/Fluids-Research/programs/P02_singularity_ladders
export OMP_NUM_THREADS=1
timeout 40000 python3 ladder_scan.py 131072 F128k_dn 900 120 0.7 -1 > ladF_128k_dn.log 2>&1 &
timeout 40000 python3 ladder_scan.py 262144 F256k_dn 900 120 0.7 -1 > ladF_256k_dn.log 2>&1 &
