#!/bin/bash
cd "$(dirname "$0")"
PARS=0 ./run.sh > run_even.log 2>&1; ./pd.sh even 0 > pd_even.log 2>&1 &
PARS=1 ./run.sh > run_odd.log 2>&1; wait; ./pd.sh odd 1 > pd_odd.log 2>&1
echo done > all.done
