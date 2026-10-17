#!/bin/bash
# Grok S66 follow-up chain (detached): waits for the L=71/50 m=1/2 runs, then m-sensitivity / reduced-N runs. NUMERICAL.
cd /workspace/GPPDiscovery5/s66; source /workspace/venv/bin/activate
while kill -0 186179 2>/dev/null || kill -0 186180 2>/dev/null; do sleep 20; done
python mori_route.py --L 7/5 --N 786 --Omega 336 --s 464/125 --m 1/8 --u0 205 --u1 68 --prec 512 --parity 0 --workers 4 --png 96 --delta_log2 1380 --out run/L75.jsonl > run/L75_m18_N786_p0.log 2>&1 &
python mori_route.py --L 71/50 --N 1006 --Omega 424 --s 763/200 --m 1/4 --u0 272 --u1 90 --prec 512 --parity 0 --workers 4 --png 96 --delta_log2 1761 --out run/L7150.jsonl --dump run/L7150_m14_N1006_p0.pkl > run/L7150_m14_N1006_p0.log 2>&1 &
wait
python mori_route.py --L 7/5 --N 800 --Omega 384 --s 464/125 --m 1/4 --u0 249 --u1 83 --prec 512 --parity 0 --workers 4 --png 96 --delta_log2 1575 --out run/L75.jsonl > run/L75_m14_N800_p0.log 2>&1 &
python mori_route.py --L 71/50 --N 1006 --Omega 424 --s 763/200 --m 1/4 --u0 272 --u1 90 --prec 512 --parity 1 --workers 4 --png 96 --delta_log2 1761 --out run/L7150.jsonl --dump run/L7150_m14_N1006_p1.pkl > run/L7150_m14_N1006_p1.log 2>&1 &
wait
echo CHAIN_DONE
