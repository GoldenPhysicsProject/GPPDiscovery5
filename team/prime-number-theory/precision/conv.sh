#!/bin/bash
# convergence series: L prec  then N list
L=$1; PREC=$2; shift 2; tag=$(echo $L | tr / _)
for N in "$@"; do ../../../../venv/bin/python legweil.py $L $N $PREC 4 60 >> runs/leg_L${tag}.jsonl; done
