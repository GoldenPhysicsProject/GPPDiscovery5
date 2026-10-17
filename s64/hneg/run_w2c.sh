#!/bin/bash
cd "$(dirname "$0")"
( for L in 1.0 1.2 1.4 1.6; do echo "PREC=400 -1 $L 64 even"; echo "PREC=400 -0.5 $L 64 even"; done
  for L in 1.8 2.0; do echo "PREC=400 -0.5 $L 64 even"; done
  for L in 1.4 1.6; do echo "PREC=400 0 $L 64 even"; done
  for L in 1.2 1.4; do echo "PREC=400 -1 $L 64 odd"; echo "PREC=400 -0.5 $L 64 odd"; done
) | xargs -P 3 -L 1 sh -c 'env $0 /workspace/venv/bin/python ht_window2.py $1 $2 $3 $4 >/dev/null 2>>w2c.err'
echo done > w2c.done
