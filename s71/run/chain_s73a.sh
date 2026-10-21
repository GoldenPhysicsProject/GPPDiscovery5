#!/bin/bash
# Grok S73 H73a: odd parity at support 2.84 (L=71/50) failed at Omega 1000 (lambda_min in (-9.3e-78,-5.9e-78), numerical bisect, CERTIFICATE FAILURE not counterexample). Each parity may use its own Omega; even certified at 1000. Retry odd at Omega 1150 (+15%). N=2.376*Omega, u0=0.607*Omega, u1=0.2*Omega, delta_log2 4900. (2026-10-21) Relaunch 19:0xZ with 3 workers: 6 workers OOM-killed at N 2732 (~3.5 GB per worker, 16 GB box).
# S66/S68: Omega 488, N 1142 gave lambda_min(T_even) ~ -4.66e-49 (NUMERICAL; S66 positivity retracted in S67).
# Scaling from S66: N ~ 2.34 Omega, u0 ~ 0.607 Omega, u1 ~ 0.2 Omega, s = 464/125. Even parity only, runs strictly sequential.
# S69 fix (07:5xZ): first launch omitted --delta_log2 (default 0) -> assert failed after ~1h assembly; now DL > -log2 d(c) (2448.05 / 2852.45), assert moved before assembly.
# Per run: certified assembly -> PSD sign test (s36 cholres at l0 = 1e-250) -> nearest-zero est -> s36 at 0.999/1.005 est -> h68b bisect if sign test fails.
cd /workspace/GPPDiscovery5/s71; source /workspace/venv/bin/activate
for P in 1; do for cfg in "1150 2732 698 230 4900"; do
  set -- $cfg; O=$1; N=$2; U0=$3; U1=$4; DL=$5; tag=L7150_m12_O${O}_p$P
  echo "$tag assembly start $(date -u +%H:%M:%S)"
  python mori_cert.py --L 71/50 --N $N --Omega $O --s 763/200 --m 1/2 --u0 $U0 --u1 $U1 --prec 512 --parity $P --workers 3 --png 96 --delta_log2 $DL --out run/cert_s73.jsonl --dump run/T_$tag.json > run/cert_$tag.log 2>&1 || { echo "$tag ASSEMBLY FAILED"; tail -5 run/cert_$tag.log; continue; }
  f=run/T_$tag.json; sha=$(sha256sum $f | cut -d' ' -f1); echo "$tag assembled sha256 $sha $(date -u +%H:%M:%S)"
  python s36_cholres.py $f 640 1 1e-250 > run/psdtest_$tag.log 2>&1; echo "$tag PSD test (T - 1e-250 I):"; grep "fac=" run/psdtest_$tag.log
  python lam_from_dump.py $f > run/lam_$tag.json; est=$(python -c "import json;print(json.load(open('run/lam_$tag.json'))['lam_min_est'])"); echo "$tag lam_est(nearest-zero, NOT lambda_min if indefinite) $est"
  if grep -q "minpiv" run/psdtest_$tag.log; then
    python s36_cholres.py $f 640 0.999,1.005 $est > run/cholres_$tag.log 2>&1; grep "fac=" run/cholres_$tag.log
    pos=$(python -c "print('%.4e'%(0.999*$est))"); ctl=$(python -c "print('%.4e'%(1.005*$est))")
    python s37_ldlt_exact.py $f $sha 640 $pos,$ctl run/ldlt_$tag.json > run/ldlt_$tag.log 2>&1; tail -3 run/ldlt_$tag.log
  else
    python h68b_bisect_vec.py $f 1 -20 -120 run/negvec_$tag.json > run/negvec_$tag.log 2>&1; tail -6 run/negvec_$tag.log
  fi
  echo "$tag done $(date -u +%H:%M:%S)"
done; done
echo "CHAIN_S73A_DONE $(date -u +%H:%M:%S)"
