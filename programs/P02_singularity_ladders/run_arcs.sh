#!/bin/bash
cd /home/user/Fluids-Research/programs/P02_singularity_ladders
timeout 3500 python3 arclength.py 2048 -1 400 > arc_down_2048.log 2>&1 &
timeout 3500 python3 arclength.py 2048 1 300 > arc_up_2048.log 2>&1 &
