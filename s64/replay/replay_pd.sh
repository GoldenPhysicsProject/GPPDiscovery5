#!/bin/bash
# Grok S64 H64a: standalone LEVEL-1 replay of the support-2.8 (L = 7/5) PD step from a stored Gram dump.
# usage: ./replay_pd.sh even|odd /path/to/atil_dump.json [outdir]
# Needs: python3 with python-flint (tested 0.9.0), numpy, mpmath. Set PY=... to override the interpreter.
# What it checks (H3 of the S35 conditional theorem), for the 384-bit midpoint matrix Atil in the dump:
#   (1) sha256 of the dump equals the value in the certificate JSON;
#   (2) s36_cholres (Arb, 640 bits, approx Cholesky + rigorous residual): lambda_min(Atil) >= l0 - ||E|| - delta_str;
#   (3) s37_ldlt_exact (exact integers, fixed-point LDL^T + exact residual + Weyl, S = 640): same, no Arb;
#   (4) NEGATIVE CONTROL: s37 at 1.02 x lam_est must FAIL (shows the test is sharp, not vacuous);
#   (5) s55_lam: Sylvester inertia (0 negative pivots) of a blocked midpoint LDL^T.
# PASS iff (1)-(5) hold and both bounds reproduce the certificate's lambda_min to 4 significant digits.
# The final window bound then follows by the certificate's bookkeeping: Q >= min(bounds) - rE - delta_str (two-block min with lamD).
set -u; P=$1; F=$2; OUT=${3:-replay_out_$P}; PY=${PY:-/workspace/venv/bin/python}
cd "$(dirname "$0")"; mkdir -p "$OUT"
case $P in even) C=cert_L140_EVEN_S61_FINAL.json;; odd) C=cert_L140_ODD_S63_FINAL.json;; *) echo "parity?"; exit 2;; esac
SHA_EXP=$($PY -c "import json;print(json.load(open('$C'))['atil_sha256'])")
LE=$($PY -c "import json;print(json.load(open('$C'))['lam_est'])")
L0=$($PY -c "print('%.4e'%($LE*0.999))"); NEG=$($PY -c "print('%.4e'%($LE*1.02))")
T0=$(date +%s)
SHA=$(sha256sum "$F" | cut -d' ' -f1)
[ "$SHA" = "$SHA_EXP" ] && echo "sha256 OK $SHA" > $OUT/sha.log || { echo "SHA MISMATCH $SHA vs $SHA_EXP" | tee $OUT/sha.log; exit 1; }
$PY s36_cholres.py "$F" 640 0.999 $LE > $OUT/cholres.log 2>&1
$PY s37_ldlt_exact.py "$F" $SHA 640 $L0 $OUT/ldlt.json > $OUT/ldlt.log 2>&1
$PY s37_ldlt_exact.py "$F" $SHA 640 $NEG $OUT/neg.json > $OUT/neg.log 2>&1
$PY s55_lam.py "$F" > $OUT/lam.log 2>&1
$PY - "$OUT" "$C" "$P" $T0 <<'PYEOF'
import sys, json, re, time
out, C, P, T0 = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]); cert = json.load(open(C))
ch = open(f'{out}/cholres.log').read(); m = re.search(r'lambda_min\(Atil\) >= ([0-9.e+-]+) OK', ch)
ld = json.load(open(f'{out}/ldlt.json')); ng = json.load(open(f'{out}/neg.json'))
lam = open(f'{out}/lam.log').read()
b1 = float(m.group(1)) if m else None; b2 = ld.get('lam_lower')
neg_fail = all(not t.get('ok', False) for t in ng['tries'])
inert = bool(re.search(r'"n_neg_pivots": 0,', lam)) and bool(re.search(r'"n_zero_pivots": 0,', lam))
ref = float(re.search(r'>= ([0-9.e+-]+)', cert['pd_impl_2_ldlt_exact_S640']).group(1))
agree = all(x is not None and abs(x / ref - 1) < 1e-3 for x in (b1, b2))
r = {'parity': P, 'cholres_lower': b1, 'ldlt_exact_lower': b2, 'neg_control_fails': neg_fail,
     'neg_control_pivot': [t.get('pivot') for t in ng['tries']], 'inertia_ok': inert, 'inertia_log': lam.strip()[-300:],
     'cert_ldlt_lower': ref, 'agree_4sig': agree, 'wall_s': time.time() - T0,
     'PASS': bool(b1 and b2 and b1 > 0 and b2 > 0 and neg_fail and agree and inert)}
json.dump(r, open(f'{out}/REPLAY_RESULT.json', 'w'), indent=1); print(json.dumps(r, indent=1))
PYEOF
echo done > $OUT/replay.done
