"""Grok S38: combine onestroke run JSON with cholres (and optional ldlt) into FINAL cert. RH OPEN."""
import json, sys
run = sys.argv[1]
ev = sys.argv[2]
od = sys.argv[3]  # may be 'NONE' if odd indefinite
out = sys.argv[4]
d = json.load(open(run))
for par, fn in (('even', ev), ('odd', od)):
    e = d.get(par)
    if e is None:
        print(par, 'MISSING in run'); continue
    if fn == 'NONE' or e.get('lam_est') is None or (e.get('lam_est') is not None and e['lam_est'] <= 0):
        e['pd_method'] = 'FAILED: Gram midpoint indefinite (lam_est<=0)'
        e['certified_lower'] = None
        print(par, 'INDEFINITE lam_est=', e.get('lam_est')); continue
    c = json.load(open(fn))
    assert c['K'] == d['K'] and 'lam_lower' in c
    e['lam0_schur'] = repr(c['lam_lower'])
    e['pd_method'] = 's38_cholres (approx Cholesky + rigorous Arb residual), prec %d' % c['prec']
    e['cholres'] = {k: c[k] for k in ('l0', 'lam_lower', 'delta_str', 'prec')}
    e['cholres']['normE'] = c['tries'][-1].get('normE')
    e['certified_lower'] = min(c['lam_lower'] - e['normE_bound'], e['lamD_lower']) - e['eps_B']
    print(par, 'certified_lower', e['certified_lower'])
d['A_cert'] = '4.4702'
d['torus'] = '../../s36/torus_L1.3_D8x2_p235_p71113_cert.json'
d['lift_note'] = 'D8x2 separating test (vs S36 D4+D20 A=4.6565)'
json.dump(d, open(out, 'w'), indent=1)
print('wrote', out)
