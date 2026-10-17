#!/bin/bash
# Grok S63 H63a: simplicity at support 2.8 (L = 7/5). Courant-Fischer restriction of the S61 even dump (K = 1497, T = 1310,
# sha256 94238400...), exactly as S62 did at 27/20; PD of M' - s2 I by two independent implementations + a failing control.
cd "$(dirname "$0")"; PY=/workspace/venv/bin/python; D=../s58/run140/atil_L140_T1310_K1497_even.json
$PY lam2_dump.py $D 768 mprime_L140_even.json > lam2_dump_L140.log 2>&1 || { echo "dump failed" > lam2.done; exit 1; }
S2=$($PY -c "import json;print(json.load(open('mprime_L140_even_meta.json'))['s2'])")
L2=$($PY -c "import json;print(json.load(open('mprime_L140_even_meta.json'))['lam2_est'])")
CTRL=$($PY -c "print('%.4e'%(float('$L2')*1.1))")
SHA=$(sha256sum mprime_L140_even.json | cut -d' ' -f1); echo "s2=$S2 lam2_est=$L2 ctrl=$CTRL sha=$SHA" > lam2_params.txt
$PY ../s36/s36_cholres.py mprime_L140_even.json 640 1.0 $S2 > pd1_lam2_L140.log 2>&1; mv mprime_L140_even_cholres_p640.json mprime_L140_even_cholres_p640_s2.json
$PY ../s36/s36_cholres.py mprime_L140_even.json 640 1.0 $CTRL > pd1_lam2_L140_ctrl.log 2>&1; mv mprime_L140_even_cholres_p640.json mprime_L140_even_cholres_p640_ctrl.json
$PY ../s37/s37_ldlt_exact.py mprime_L140_even.json $SHA 640 $S2,$CTRL mprime_L140_even_ldltexact_S640.json > pd2_lam2_L140.log 2>&1
echo done > lam2.done
