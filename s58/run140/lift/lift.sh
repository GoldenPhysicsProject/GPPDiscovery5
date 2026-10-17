#!/bin/bash
# Grok S58 H58a step 1: BUFFERED lift at L = 7/5 (support 2.8; comb now includes n = 16).
# Bundle LP (s55_bundle_lp, EPS = 0.0056: every character |k.logp| >= 2L + 0.0056), warm-started from the S55 L = 1.38 bundles,
# then two B&B implementations per part (float s34_bb_D4 SUBSET + libm-free s35_interval_bb), separable merge + group-locality check.
set -e; cd "$(dirname "$0")"; PY=/workspace/venv/bin/python; R=/workspace/GPPDiscovery5
export EPS=0.0056
run_part() { P=$1; D=$2; PR=$3; INIT=$4; MB=$5; MG=$6; BUD=$7
  mkdir -p bb_$P; ( cd bb_$P
    $PY $R/s55/s55_bundle_lp.py 1.4 $D $BUD L140_D${D}_${P}_buf $PR $INIT > lp.log 2>&1
    $PY $R/s54/mk_src.py bundle_L140_D${D}_${P}_buf.json src.json > mksrc.log 2>&1
    SRC=src.json SUBSET=1 MAXBOXES=$MB $PY $R/s34/s34_bb_D4.py L1.4_D$D $MG > bb.log 2>&1
    TG=$($PY -c "import json;print(json.load(open('torus_L1.4_D${D}_cert.json'))['A_cert'])")
    $PY $R/s35/s35_interval_bb.py torus_L1.4_D${D}_cert.json $TG 20 7/5 > iv.log 2>&1 ) ; }
run_part p235 9 2,3,5 $R/s55/lift/bundle_L138_D9_p235_buf.json 12000000 0.003 1100 &
run_part p71113 8 7,11,13 $R/s55/lift/bundle_L138_D8_p71113_buf.json 40000000 0.004 1000 &
wait
grep -q '"ok": true' bb_p235/iv.log && grep -q '"ok": true' bb_p71113/iv.log || { echo "INTERVAL B&B FAILED"; exit 3; }
$PY $R/s34/s34_merge_torus2.py bb_p235/torus_L1.4_D9_cert.json bb_p71113/torus_L1.4_D8_cert.json 7/5 torus_L1.40_buf_D9D8_cert.json
$PY $R/s35/s35_check_merge.py torus_L1.40_buf_D9D8_cert.json bb_p235/torus_L1.4_D9_cert.json bb_p71113/torus_L1.4_D8_cert.json:3 > check_merge.log 2>&1
cat check_merge.log
echo LIFT DONE
