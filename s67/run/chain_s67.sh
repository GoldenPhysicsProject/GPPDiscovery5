#!/bin/bash
# Grok S67 detached chain: certified-assembly reruns of L=71/50 m=1/2 N=1272 (both parities in parallel), delta=2^-2221 <= d(c).
cd /workspace/GPPDiscovery5/s67; source /workspace/venv/bin/activate
for p in 0 1; do
python mori_cert.py --L 71/50 --N 1272 --Omega 536 --s 763/200 --m 1/2 --u0 308 --u1 102 --prec 512 --parity $p --workers 8 --png 96 --delta_log2 2221 --out run/cert_L7150.jsonl --dump run/T_L7150_m12_N1272_p$p.json > run/cert_L7150_p$p.log 2>&1 &
done
wait; echo "ASSEMBLY_DONE $(date -u +%H:%M:%S)"
