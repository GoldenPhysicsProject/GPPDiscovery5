#!/bin/bash
# Grok S67 detached PD chain: waits for each certified-assembly dump, then two independent PD certificates
# (s36_cholres: Arb residual; s37_ldlt_exact: exact-integer LDL^T) at l0 just below lambda_est, plus negative controls above it.
cd /workspace/GPPDiscovery5/s67; source /workspace/venv/bin/activate
declare -A POS=( [0]=1.850e-81 [1]=2.858e-77 ); declare -A CTL=( [0]=1.86e-81 [1]=2.87e-77 ); declare -A EST=( [0]=1.85047476363e-81 [1]=2.85827975703e-77 )
for p in 0 1; do
  f=run/T_L7150_m12_N1272_p$p.json
  while ! grep -q "\"parity\": $p" run/cert_L7150.jsonl 2>/dev/null; do sleep 15; done
  sha=$(sha256sum $f | cut -d' ' -f1); echo "p$p sha256 $sha $(date -u +%H:%M:%S)"
  fac=$(python -c "print(${POS[$p]}/${EST[$p]})"); facc=$(python -c "print(${CTL[$p]}/${EST[$p]})")
  python s36_cholres.py $f 640 $fac ${EST[$p]} > run/cholres_p$p.log 2>&1
  cp run/T_L7150_m12_N1272_p${p}_cholres_p640.json run/cholres_p${p}_pos.json 2>/dev/null
  python s36_cholres.py $f 640 $facc ${EST[$p]} > run/cholres_p${p}_ctrl.log 2>&1
  cp run/T_L7150_m12_N1272_p${p}_cholres_p640.json run/cholres_p${p}_ctrl.json 2>/dev/null
  python s37_ldlt_exact.py $f $sha 640 ${POS[$p]},${CTL[$p]} run/ldlt_p$p.json > run/ldlt_p$p.log 2>&1 &
done
wait; echo PD_CHAIN_DONE $(date -u +%H:%M:%S)
