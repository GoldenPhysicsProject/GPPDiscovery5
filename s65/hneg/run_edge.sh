#!/bin/bash
# Grok S65 H65c edge chain: refine t_c of the gamma_34/35 pair; full scans + complex zeros + K=64 windows at t=-0.4,-0.45
cd "$(dirname "$0")"
PY=/workspace/venv/bin/python
for t in -0.37 -0.38 -0.385 -0.39; do echo $t; done | xargs -P 3 -I{} sh -c "$PY tc_track.py {} >> tc.jsonl 2>>tc.err"
for t in -0.4 -0.45; do echo $t; done | xargs -P 2 -I{} sh -c "$PY ht_scan.py {} 300 0.02 scan_t{}.json > scan_t{}.log 2>&1 && $PY ht_complex.py {} > cplx_t{}.log 2>&1"
(for t in -0.4 -0.45; do for L in 1.2 1.3 1.4 1.6; do for p in even odd; do echo "$t $L 64 $p"; done; done; done) | xargs -P 3 -n 4 sh -c "$PY ht_window2.py \$0 \$1 \$2 \$3 >/dev/null 2>>w.err"
echo done > edge.done
