#!/bin/bash
# Grok S56 -> S57: L = 69/50 ODD parity, T = 1160, K = 1305 (beta = 0.2533), CERTIFIED buffered lift A <= 4.9641.
cd "$(dirname "$0")"
export T=1160 K=1305
TORUS=torus_L1.38_buf_D9D8_cert.json PARS=1 W=4 ./run.sh > run_odd.log 2>&1 && ./pd.sh odd 1 > pd_odd.log 2>&1
/workspace/venv/bin/python ../../s55/s55_lam.py atil_L138_T1160_K1305_odd.json > lam_odd.log 2>&1
echo done > odd.done
