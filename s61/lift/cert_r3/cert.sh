#!/bin/bash
# Grok S61 H61b: certify the round-3 D10 {2,3,5} bundle (LP sample 3.151787) exactly as S58: float B&B (s34, SUBSET) + libm-free interval B&B (s35),
# then merge with the S58 certified {7,11,13} D8 part and check group-locality.
set -e; cd "$(dirname "$0")"; PY=/workspace/venv/bin/python; R=/workspace/GPPDiscovery5; export EPS=0.0056 OMP_NUM_THREADS=1
$PY $R/s54/mk_src.py bundle_L140_D10_p235_buf.json src.json > mksrc.log 2>&1
SRC=src.json SUBSET=1 MAXBOXES=12000000 nice -n 3 $PY $R/s34/s34_bb_D4.py L1.4_D10 0.003 > bb.log 2>&1
TG=$($PY -c "import json;print(json.load(open('torus_L1.4_D10_cert.json'))['A_cert'])")
$PY $R/s35/s35_interval_bb.py torus_L1.4_D10_cert.json $TG 20 7/5 > iv.log 2>&1
grep -q '"ok": true' iv.log || { echo "INTERVAL B&B FAILED" > cert.done; exit 3; }
P71=$R/s58/run140/lift/bb_p71113/torus_L1.4_D8_cert.json
$PY $R/s34/s34_merge_torus2.py torus_L1.4_D10_cert.json $P71 7/5 torus_L1.40_buf_D10D8_cert.json > merge.log 2>&1
$PY $R/s35/s35_check_merge.py torus_L1.40_buf_D10D8_cert.json torus_L1.4_D10_cert.json $P71:3 > check_merge.log 2>&1
cat check_merge.log; echo done > cert.done
