#!/bin/bash
# Grok S67-resumed detached chain (2026-10-21): fixed mori_cert.py (midpoint+derivative radii for Legendre and
# spherical-Bessel recurrences), L=71/50 m=1/2 N=1272, parities STRICTLY SEQUENTIAL (15 GB box, shared).
# Per parity: certified assembly -> NUMERICAL lam estimate -> s36 cholres at l0 = 0.999*est (positive) and 1.005*est (control, must fail)
# -> s37 exact-integer LDL^T at the same two l0.
cd /workspace/GPPDiscovery5/s67; source /workspace/venv/bin/activate
for p in 0 1; do
  echo "p$p assembly start $(date -u +%H:%M:%S)"
  python mori_cert.py --L 71/50 --N 1272 --Omega 536 --s 763/200 --m 1/2 --u0 308 --u1 102 --prec 512 --parity $p --workers 8 --png 96 --delta_log2 2221 --out run/cert_L7150_r.jsonl --dump run/T_L7150_m12_N1272_p$p.json > run/cert_L7150_r_p$p.log 2>&1 || { echo "p$p ASSEMBLY FAILED"; exit 1; }
  f=run/T_L7150_m12_N1272_p$p.json; sha=$(sha256sum $f | cut -d' ' -f1); echo "p$p assembled sha256 $sha $(date -u +%H:%M:%S)"
  python lam_from_dump.py $f > run/lam_p$p.json; est=$(python -c "import json;print(json.load(open('run/lam_p$p.json'))['lam_min_est'])"); echo "p$p lam_est $est"
  python s36_cholres.py $f 640 0.999,1.005 $est > run/cholres_r_p$p.log 2>&1; echo "p$p cholres done $(date -u +%H:%M:%S)"; grep "fac=" run/cholres_r_p$p.log
  pos=$(python -c "print('%.4e'%(0.999*$est))"); ctl=$(python -c "print('%.4e'%(1.005*$est))")
  python s37_ldlt_exact.py $f $sha 640 $pos,$ctl run/ldlt_r_p$p.json > run/ldlt_r_p$p.log 2>&1; echo "p$p ldlt done $(date -u +%H:%M:%S)"; tail -4 run/ldlt_r_p$p.log
done
echo "CHAIN_S67R_DONE $(date -u +%H:%M:%S)"
