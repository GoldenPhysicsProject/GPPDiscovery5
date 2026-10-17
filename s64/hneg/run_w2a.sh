#!/bin/bash
cd "$(dirname "$0")"
(for t in -0.5 -1; do for L in 1.4 1.5 1.6 1.7 1.8 2.0; do for p in even odd; do echo "$t $L 32 $p"; done; done; done; for L in 1.6 2.0; do echo "0 $L 32 even"; done) | xargs -P 4 -n 4 sh -c '/workspace/venv/bin/python ht_window2.py $0 $1 $2 $3 >/dev/null 2>>w2a.err'
echo done > w2a.done
