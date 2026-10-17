#!/bin/bash
# Grok S60: re-run of S58 even L=7/5 T=1310 K=1497 with SPHJ_BOOST=2.0. The S58 run (BOOST 1.75) was VOID:
# sph_j_all backward recurrence at x = T*L = 1834, D = 2994 lost all precision (j_0 ball +/- 2e731, rE = inf; see void_boost175/).
cd "$(dirname "$0")"; PY=/workspace/venv/bin/python
export SPHJ_BOOST=2.0 T=1310 K=1497
PARS=0 W=${W:-3} ./run.sh > run_even.log 2>&1 || { echo "run failed" > even.done; exit 1; }
grep -q 'rE=inf' run_even.log && { echo "rE=inf again" > even.done; exit 1; }
./pd.sh even 0 > pd_even.log 2>&1
$PY ../../s55/s55_lam.py atil_L140_T${T}_K${K}_even.json > lam_even.log 2>&1
echo done > even.done
