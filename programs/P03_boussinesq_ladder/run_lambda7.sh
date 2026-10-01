#!/bin/bash
# matched λ7 scans and the resolved endpoint run; every job resumes from its saved states after a restart
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /home/user/Fluids-Research/programs/P03_boussinesq_ladder
SAVE_EVERY=1 NK_TOL=7e-13 SS_IN=-16 nohup python3 ipm_scan.py ipm_scan_s16_z7.2260.npy 0.1383891503 7.196 7.4961 0.03 -20 p7 0.0125 > /dev/null 2>&1 &
SAVE_EVERY=1 NK_TOL=1.5e-12 HS_IN=0.00625 nohup python3 ipm_scan.py ipm_deep_e2_lam0.1376273.npy 0.1376273 7.196 7.4961 0.06 -20 h7 0.00625 > /dev/null 2>&1 &
SAVE_EVERY=1 NK_TOL=8e-13 NB=48 SS_IN=-16 nohup python3 ipm_scan.py ipm_scan_s16_z7.2260.npy 0.1383891503 7.196 7.4961 0.06 -20 n7 0.0125 > /dev/null 2>&1 &
IPM_SSTART=-20 IPM_SSTART_IN=-20 nohup python3 ipm_deep.py 0.1376273 0.07 0.25 e3 ipm_deep_e2_lam0.1376273.npy 32 0.00625 32 0.00625 1e-10 > /dev/null 2>&1 &
