#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
until ! pgrep -f "wkb_sens.py" > /dev/null; do sleep 30; done
python3 wkb_branch.py 3 > wkb_branch_D3.out 2>&1
python3 wkb_branch.py 2 > wkb_branch_D2.out 2>&1
