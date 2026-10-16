"""Grok S36: combine the schurboost run JSON (Gram assembly, rE, tails) with the cholres PD certificates
(s36_cholres.py, replaces the failed ball-Schur step) into one certificate JSON; certified_lower recomputed as
min(lam_chol - rE, lamD) - eps_B with lam_chol = l0 - ||E|| - delta_str (already subtracted). RH OPEN."""
import json, sys
run, ev, od, out = sys.argv[1:5]
d = json.load(open(run))
for par, fn in (('even', ev), ('odd', od)):
    c = json.load(open(fn)); e = d[par]
    assert c['K'] == d['K'] and 'lam_lower' in c
    e['lam0_schur'] = repr(c['lam_lower']); e['pd_method'] = 's36_cholres (approx Cholesky + rigorous Arb residual), prec %d' % c['prec']
    e['cholres'] = {k: c[k] for k in ('l0', 'lam_lower', 'delta_str', 'prec')}; e['cholres']['normE'] = c['tries'][-1]['normE']
    e['certified_lower'] = min(c['lam_lower'] - e['normE_bound'], e['lamD_lower']) - e['eps_B']
    print(par, 'certified_lower', e['certified_lower'])
json.dump(d, open(out, 'w'), indent=1)
