#!/bin/bash
# Grok S61: ODD parity of support 2.8 (L = 7/5). Mirrors chain_even_r2.sh (SPHJ_BOOST = 2.0, W = 3), but uses the S61
# B&B-certified {2,3,5} D10 + S58 {7,11,13} D8 lift (TORUS env), T chosen from A_cert for beta_odd target 0.25 (S57 odd used 0.2533).
cd "$(dirname "$0")"; PY=/workspace/venv/bin/python
export SPHJ_BOOST=2.0 TORUS=${TORUS:?}
read T K <<< $($PY -c "
import json,math
A=float(json.load(open('$TORUS'))['A_cert']); b=float('${BETA:-0.25}')
T=math.ceil(2*math.pi*math.exp(A+b)+1); K=math.ceil(0.8152*1.4*T)+1
print(T,K)")
export T K
echo "TORUS=$TORUS A_cert=$($PY -c "import json;print(json.load(open('$TORUS'))['A_cert'])") T=$T K=$K beta_target=${BETA:-0.25}" > params_odd.txt
PARS=1 W=${W:-3} ./run.sh > run_odd.log 2>&1 || { echo "run failed" > odd.done; exit 1; }
grep -q 'rE=inf' run_odd.log && { echo "rE=inf" > odd.done; exit 1; }
./pd.sh odd 1 > pd_odd.log 2>&1
$PY ../../s55/s55_lam.py atil_L140_T${T}_K${K}_odd.json > lam_odd.log 2>&1
echo done > odd.done
