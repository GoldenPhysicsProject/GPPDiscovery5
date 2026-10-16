#!/bin/bash
# Run cholres + ldlt on dumps after assembly. Skip parity if lam_est<=0.
set -euo pipefail
cd "$(dirname "$0")"
PY=/workspace/venv/bin/python
CERT=cert_L13_T700_K740_ng78_D8x2_p384.json
PREFIX=atil_L13_T700_K740_D8x2
echo "=== PD start $(date -u +%Y-%m-%dT%H:%M:%SZ) ===" | tee -a pd.log

sha256sum ${PREFIX}_even.json ${PREFIX}_odd.json > atil_sha256.txt
cat atil_sha256.txt | tee -a pd.log

EVEN_LAM=$($PY -c "import json; d=json.load(open('$CERT')); print(d['even']['lam_est'])")
ODD_LAM=$($PY -c "import json; d=json.load(open('$CERT')); print(d['odd']['lam_est'])")
echo "even lam_est=$EVEN_LAM odd lam_est=$ODD_LAM" | tee -a pd.log

EV_CHOL=NONE
OD_CHOL=NONE

if $PY -c "import sys; sys.exit(0 if float('$EVEN_LAM')>0 else 1)"; then
  echo "=== cholres EVEN ===" | tee -a pd.log
  $PY ../../s38/s38_cholres.py ${PREFIX}_even.json 640 0.999,0.99,0.95 "$EVEN_LAM" 64 2>&1 | tee cholres_even_p640.log
  EV_CHOL=${PREFIX}_even_cholres_p640.json
  SHA=$(awk '/even.json/{print $1}' atil_sha256.txt)
  L0=$($PY -c "import json; print('%.6g'%json.load(open('$EV_CHOL'))['l0'])")
  echo "=== ldlt EVEN l0=$L0 ===" | tee -a pd.log
  $PY ../../s38/s38_ldlt_exact.py ${PREFIX}_even.json "$SHA" 640 "$L0" 2>&1 | tee ldlt_even.log || true
else
  echo "EVEN indefinite — skip PD" | tee -a pd.log
fi

if $PY -c "import sys; sys.exit(0 if float('$ODD_LAM')>0 else 1)"; then
  echo "=== cholres ODD ===" | tee -a pd.log
  $PY ../../s38/s38_cholres.py ${PREFIX}_odd.json 640 0.999,0.99,0.95 "$ODD_LAM" 64 2>&1 | tee cholres_odd_p640.log
  OD_CHOL=${PREFIX}_odd_cholres_p640.json
  SHA=$(awk '/odd.json/{print $1}' atil_sha256.txt)
  L0=$($PY -c "import json; print('%.6g'%json.load(open('$OD_CHOL'))['l0'])")
  echo "=== ldlt ODD l0=$L0 ===" | tee -a pd.log
  $PY ../../s38/s38_ldlt_exact.py ${PREFIX}_odd.json "$SHA" 640 "$L0" 2>&1 | tee ldlt_odd.log || true
else
  echo "ODD indefinite — skip PD (CONJECTURE CONFIRMED: lift-induced)" | tee -a pd.log
fi

$PY ../../s38/s38_finalize.py "$CERT" "$EV_CHOL" "$OD_CHOL" cert_L13_T700_K740_D8x2_FINAL.json 2>&1 | tee -a pd.log
echo "=== PD done $(date -u +%Y-%m-%dT%H:%M:%SZ) ===" | tee -a pd.log
