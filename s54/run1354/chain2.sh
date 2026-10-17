#!/bin/bash
# wait for chain.sh, then PD checks (run after each parity: even PD overlaps odd assembly only if memory allows; keep sequential)
cd "$(dirname "$0")"
while [ ! -f cert_L1354_T800_K884_ng92_p384_par0.json ]; do sleep 10; done
./pd.sh even 0 > pd_even.log 2>&1
while [ ! -f chain.done ]; do sleep 10; done
./pd.sh odd 1 > pd_odd.log 2>&1
echo ALL DONE > all.done
