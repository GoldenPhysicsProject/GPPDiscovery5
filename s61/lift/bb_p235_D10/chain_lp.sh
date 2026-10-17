#!/bin/bash
# Grok S61 H61b: D10 {2,3,5} buffered bundle LP at L=7/5 in memory-capped rounds (each round restarts the cutting-plane
# point set, keeping RSS < ~2.2 GB so the S58 certification run is never starved). Warm start = previous round's bundle.
cd "$(dirname "$0")"; export EPS=0.0056 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/workspace/venv/bin/python; R=/workspace/GPPDiscovery5
for k in 2 3 4 5; do
  cp bundle_L140_D10_p235_buf.json init_r$k.json
  nice -n 5 $PY $R/s55/s55_bundle_lp.py 1.4 10 ${BUD:-560} L140_D10_p235_buf 2,3,5 init_r$k.json > lp_round$k.log 2>&1
  [ -f stop ] && break
done
echo done > lp.done
