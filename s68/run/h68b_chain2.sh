#!/bin/bash
cd /workspace/GPPDiscovery5/s68; source /workspace/venv/bin/activate
python h68b_bisect_vec.py ../s67/run/T_L7150_m12_N1272_p1.json 1 -20 -40 run/h68b_vec_p1.json > run/h68b_vec_p1.log 2>&1
python h68b_bisect_vec.py ../s67/run/s66check/L75_m12_N1142_p0.json 0 -30 -60 run/h68b_vec_L75_m12_p0.json > run/h68b_vec_L75_m12_p0.log 2>&1
python h68b_bisect_vec.py ../s67/run/s66check/L75_m14_N898_p0.json 0 -30 -60 run/h68b_vec_L75_m14_p0.json > run/h68b_vec_L75_m14_p0.log 2>&1
echo H68B2_DONE $(date -u +%H:%M:%S)
