#!/bin/bash
# Grok S68 H68a detached chain (2026-10-21): Mori VIII log 15 replay with fixed S67 mori_cert.
# Mori VIII parameters: L = log(15)/2, N = 1664 (832/parity), Omega = 736, s = 33/8, m = 1/2, u0 = 428, u1 = 142, delta = 2^-2901.
# Parities STRICTLY SEQUENTIAL (shared box, ~4 GB free). Per parity: certified assembly -> NUMERICAL nearest-zero estimate
# -> s36 cholres at 0.999*est (must pass for PSD) and 1.005*est (control, must fail) -> s37 exact LDL^T same l0.
cd /workspace/GPPDiscovery5/s68; source /workspace/venv/bin/activate
for p in 0 1; do
  echo "p$p assembly start $(date -u +%H:%M:%S)"
  python mori_cert.py --L log15/2 --N 1664 --Omega 736 --s 33/8 --m 1/2 --u0 428 --u1 142 --prec 512 --parity $p --workers 6 --png 96 --delta_log2 2901 --out run/cert_log15.jsonl --dump run/T_log15_p$p.json > run/cert_log15_p$p.log 2>&1 || { echo "p$p ASSEMBLY FAILED"; tail -5 run/cert_log15_p$p.log; continue; }
  true
  f=run/T_log15_p$p.json; sha=$(sha256sum $f | cut -d' ' -f1); echo "p$p assembled sha256 $sha $(date -u +%H:%M:%S)"
  python lam_from_dump.py $f > run/lam_log15_p$p.json; est=$(python -c "import json;print(json.load(open('run/lam_log15_p$p.json'))['lam_min_est'])"); echo "p$p lam_est(nearest-zero) $est $(date -u +%H:%M:%S)"
  python s36_cholres.py $f 640 0.999,1.005 $est > run/cholres_log15_p$p.log 2>&1; echo "p$p cholres done $(date -u +%H:%M:%S)"; grep "fac=" run/cholres_log15_p$p.log
  pos=$(python -c "print('%.4e'%(0.999*$est))"); ctl=$(python -c "print('%.4e'%(1.005*$est))")
  python s37_ldlt_exact.py $f $sha 640 $pos,$ctl run/ldlt_log15_p$p.json > run/ldlt_log15_p$p.log 2>&1; echo "p$p ldlt done $(date -u +%H:%M:%S)"; tail -3 run/ldlt_log15_p$p.log
done
echo "CHAIN_S68A_DONE $(date -u +%H:%M:%S)"
