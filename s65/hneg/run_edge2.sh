#!/bin/bash
# Grok S65 H65c edge chain 2: near-t_c behaviour (t = -0.39, -0.3875) and lower L at t = -0.4
cd "$(dirname "$0")"
PY=/workspace/venv/bin/python
while [ ! -f edge.done ]; do sleep 20; done
for t in -0.39 -0.3875; do echo $t; done | xargs -P 2 -I{} sh -c "$PY ht_scan.py {} 300 0.02 scan_t{}.json > scan_t{}.log 2>&1 && $PY ht_complex.py {} > cplx_t{}.log 2>&1"
(echo "-0.4 1.0 64 even"; echo "-0.4 1.1 64 even"; echo "-0.4 1.1 64 odd"; for t in -0.39 -0.3875; do for L in 1.2 1.4; do for p in even odd; do echo "$t $L 64 $p"; done; done; done) | xargs -P 3 -n 4 sh -c "$PY ht_window2.py \$0 \$1 \$2 \$3 >/dev/null 2>>w.err"
echo done > edge2.done
