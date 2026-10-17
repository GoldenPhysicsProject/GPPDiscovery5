#!/bin/bash
# Grok S58: detached chain. lift (LP + 2x B&B) -> choose T from A_cert (beta_even target 0.15) -> even assembly -> PD (2 impls + control) -> inertia.
cd "$(dirname "$0")"; PY=/workspace/venv/bin/python
( cd lift && ./lift.sh ) > lift.log 2>&1 || { echo "lift failed" > even.done; exit 1; }
cp lift/torus_L1.40_buf_D9D8_cert.json .
read T K <<< $($PY -c "
import json,math
A=float(json.load(open('torus_L1.40_buf_D9D8_cert.json'))['A_cert']); b=float('${BETA:-0.15}')
T=math.ceil(2*math.pi*math.exp(A+b)+1); K=math.ceil(0.8152*1.4*T)+1
print(T,K)")
echo "A_cert=$(grep -o '\"A_cert\": *\"[0-9.]*\"' torus_L1.40_buf_D9D8_cert.json) T=$T K=$K" > params.txt
export T K
PARS=0 W=${W:-3} ./run.sh > run_even.log 2>&1 && ./pd.sh even 0 > pd_even.log 2>&1
$PY ../../s55/s55_lam.py atil_L140_T${T}_K${K}_even.json > lam_even.log 2>&1
echo done > even.done
