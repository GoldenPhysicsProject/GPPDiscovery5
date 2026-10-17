#!/bin/bash
cd "$(dirname "$0")"
PARS=0 ./run.sh > run_even.log 2>&1 && ./pd.sh even 0 > pd_even.log 2>&1
/workspace/venv/bin/python ../../s55/s55_lam.py atil_L138_T1020_K1148_even.json > lam_even.log 2>&1
echo done > even.done
