"""Grok S65 H65c: track the real-zero pair of H_t near x = 2*111.03, 2*111.87 (zeta gamma_34, gamma_35) as t decreases;
count sign changes on [220.5, 225] with step 0.002 (NUMERICAL, same non-rigorous quadrature as S64 ht_scan.py).
usage: tc_track.py t  -> prints t, #sign changes, min |H| location"""
import sys, json
tval = sys.argv[1]
sys.argv = ['x', tval, '300', '0.02', '/dev/null']
src = open('ht_scan.py').read().split('t0 = time.time()')[0]; exec(src)
xs = [220.5 + 0.002*i for i in range(int(4.5/0.002)+1)]
vals = [Ht(repr(x)) for x in xs]
sg = [1 if v > 0 else -1 for v in vals]
ch = [round(xs[i], 3) for i in range(len(xs)-1) if sg[i] != sg[i+1]]
lv = [abs(float(v.mid())) for v in vals]
print(json.dumps({"t": tval, "sign_changes_x": ch, "n": len(ch)}), flush=True)
