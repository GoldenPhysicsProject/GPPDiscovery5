#!/bin/bash
# Grok S56 H56a step 1: B&B-certify the S55 BUFFERED L=69/50 lift ({2,3,5} D9 + {7,11,13} D8, eps=0.0056).
# Two B&B implementations per part: float s34_bb_D4 (SUBSET mode) and libm-free s35_interval_bb; separable merge + group check.
set -e; cd "$(dirname "$0")"; PY=/workspace/venv/bin/python; R=../..
run_part() { P=$1; KEY=$2; MB=$3; MG=$4
  mkdir -p bb_$P; ( cd bb_$P
    $PY $R/../s54/mk_src.py ../bundle_L138_${KEY#L1.38_}_${P}_buf.json src.json
    SRC=src.json SUBSET=1 MAXBOXES=$MB $PY $R/../s34/s34_bb_D4.py $KEY $MG > bb.log 2>&1
    TG=$($PY -c "import json;print(json.load(open('torus_${KEY}_cert.json'))['A_cert'])")
    $PY $R/../s35/s35_interval_bb.py torus_${KEY}_cert.json $TG 20 69/50 > iv.log 2>&1 ) ; }
run_part p235 L1.38_D9 12000000 0.003 &
run_part p71113 L1.38_D8 40000000 0.004 &
wait
$PY $R/s34/s34_merge_torus2.py bb_p235/torus_L1.38_D9_cert.json bb_p71113/torus_L1.38_D8_cert.json 69/50 torus_L1.38_buf_D9D8_cert.json
$PY $R/s35/s35_check_merge.py torus_L1.38_buf_D9D8_cert.json bb_p235/torus_L1.38_D9_cert.json bb_p71113/torus_L1.38_D8_cert.json:3 > check_merge.log 2>&1
echo LIFT DONE
