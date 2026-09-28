#!/bin/bash
cd /home/user/Fluids-Research/programs/P02_singularity_ladders
export OMP_NUM_THREADS=1
timeout 30000 python3 ladder_scan.py 16384 F16k_dn 2500 120 0.7 -1 > ladF_16k_dn.log 2>&1 &
timeout 30000 python3 ladder_scan.py 16384 F16k_up 2500 120 0.7 1 > ladF_16k_up.log 2>&1 &
timeout 30000 python3 ladder_scan.py 65536 F64k_dn 600 120 0.7 -1 > ladF_64k_dn.log 2>&1 &
