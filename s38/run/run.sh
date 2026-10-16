#!/bin/bash
# Grok S38 H38a: separating test — L=1.3 with D8x2 lift A<=4.4702 (vs S36 D4+D20 A<=4.6565).
# Question: does odd Gram go indefinite under the tighter D8x2 lift?
set -euo pipefail
cd "$(dirname "$0")"
T=${T:-700}; K=${K:-740}; NG=${NG:-78}; W=${W:-2}
echo "=== H38a start $(date -u +%Y-%m-%dT%H:%M:%SZ) T=$T K=$K ng=$NG workers=$W ===" | tee -a run.log
/workspace/venv/bin/python ../../s38/s38_onestroke_schurboost.py \
  --L 13/10 --T "$T" --K "$K" --ng "$NG" --w 1/2 --b 2/5 --prec 384 \
  --workers "$W" \
  --torus ../../s36/torus_L1.3_D8x2_p235_p71113_cert.json \
  --schur_bs '' \
  --dump_atil atil_L13_T${T}_K${K}_D8x2 \
  --out cert_L13_T${T}_K${K}_ng${NG}_D8x2_p384.json \
  2>&1 | tee -a run.log
echo "=== H38a assembly done $(date -u +%Y-%m-%dT%H:%M:%SZ) ===" | tee -a run.log
