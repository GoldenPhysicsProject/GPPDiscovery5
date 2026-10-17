#!/bin/bash
# Grok S54: PD step for one parity dump (two independent implementations + negative control that must FAIL).
# usage: ./pd.sh even|odd PARNUM
set -e; cd "$(dirname "$0")"; PY=/workspace/venv/bin/python; P=$1; N=$2
F=atil_L1354_T800_K884_$P.json; C=cert_L1354_T800_K884_ng92_p384_par$N.json
LE=$($PY -c "import json;print(json.load(open('$C'))['$P']['lam_est'])")
SHA=$(sha256sum $F | cut -d' ' -f1); echo "$SHA  $F" >> atil_sha256.txt
L0=$($PY -c "print('%.4e'%($LE*0.999))"); NEG=$($PY -c "print('%.4e'%($LE*1.02))")
echo "lam_est=$LE l0=$L0 neg=$NEG"
$PY ../../s36/s36_cholres.py $F 640 0.999 $LE > cholres_${P}_p640.log 2>&1 || true
$PY ../../s37/s37_ldlt_exact.py $F $SHA 640 $L0 ${P}_ldltexact_S640.json > ldlt_${P}.log 2>&1 || true
$PY ../../s37/s37_ldlt_exact.py $F $SHA 640 $NEG neg_${P}_ldltexact_S640.json > neg_ldlt_${P}.log 2>&1 || true
echo PD $P DONE
