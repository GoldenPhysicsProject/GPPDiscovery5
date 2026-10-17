#!/bin/bash
cd "$(dirname "$0")"
( for ps in 0.05 -0.05 0.2; do echo "PSHIFT=$ps PREC=256 -1 1.5 32 even"; done
  for t in -1.5 -2 -3 -4; do for L in 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do echo "PSHIFT=0 PREC=256 $t $L 24 even"; echo "PSHIFT=0 PREC=256 $t $L 24 odd"; done; done
  for L in 1.6 1.7 1.8; do echo "PSHIFT=0 PREC=400 -0.5 $L 64 even"; echo "PSHIFT=0 PREC=400 -1 $L 64 even"; echo "PSHIFT=0 PREC=400 0 $L 64 even"; done
) | xargs -P 6 -L 1 sh -c 'env $0 $1 /workspace/venv/bin/python ht_window2.py $2 $3 $4 $5 >/dev/null 2>>w2b.err'
echo done > w2b.done
