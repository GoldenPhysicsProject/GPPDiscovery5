#!/bin/bash
# even then odd, sequential
cd "$(dirname "$0")"
PARS=0 ./run.sh > run_even.log 2>&1
PARS=1 ./run.sh > run_odd.log 2>&1
echo CHAIN DONE > chain.done
