#!/bin/bash
# Grok S54 H54a step 1: certify the L=69/50 D8x2 torus lift. Two B&B implementations per part:
# float s34_bb_D4 (SUBSET mode) and the libm-free s35_interval_bb. Then separable merge + group-locality check.
set -e; cd "$(dirname "$0")"; PY=/workspace/venv/bin/python; B=../../
for P in p235 p71113; do
  mkdir -p bb_$P; ( cd bb_$P
    $PY ../../mk_src.py ../bundle_L138_D8_$P.json src_L138_D8_$P.json
    MB=12000000; MG=0.003; [ $P = p71113 ] && MB=40000000 && MG=0.004
    SRC=src_L138_D8_$P.json SUBSET=1 MAXBOXES=$MB $PY $B/../s34/s34_bb_D4.py L1.38_D8 $MG > bb.log 2>&1
    TG=$($PY -c "import json;print(json.load(open('torus_L1.38_D8_cert.json'))['A_cert'])")
    $PY $B/../s35/s35_interval_bb.py torus_L1.38_D8_cert.json $TG 20 69/50 > iv.log 2>&1 ) &
done
wait
$PY $B/s34/s34_merge_torus2.py bb_p235/torus_L1.38_D8_cert.json bb_p71113/torus_L1.38_D8_cert.json 69/50 torus_L1.38_D8x2_p235_p71113_cert.json
$PY $B/s35/s35_check_merge.py torus_L1.38_D8x2_p235_p71113_cert.json bb_p235/torus_L1.38_D8_cert.json bb_p71113/torus_L1.38_D8_cert.json:3 > check_merge.log 2>&1
echo LIFT DONE
