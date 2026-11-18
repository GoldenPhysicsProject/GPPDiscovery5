"""Round 14E fits: zeta lambda_min (W and S^perp) at small L where it is above the float floor; fakes' lambda_min on own S^perp vs L."""
import json, math, numpy as np
from gap import *
out = dict(zeta=[], fakes=[])
Z = md.zeta()
for F in [45, 90]:
    for L in [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75]:
        M = int(math.ceil(F*L/math.pi))+8; r = analyse('ZETA', Z, L, M, F, 1e-8, W=wmatrix(Z, L, M), widths=(1.5, 1.0), kf=1.1)
        out['zeta'].append({k: r[k] for k in ['F', 'L', 'M', 'lam_min_W', 'lam_min_perp', 'lam2_perp', 'dim_perp', 'normW']}); print(out['zeta'][-1], flush=True)
for name, mo, Ls, F, kw in [('EPS', md.eps(), [3.5, 4.0, 5.0, 6.0], 60, {}), ('FTH', md.fth(0.2), [2.0, 3.0, 4.0, 5.0], 45, {}),
                            ('DH', md.dh(), [3.0, 4.0, 5.0, 6.0], 125, dict(widths=(1.5, 1.0), kf=1.0))]:
    for L in Ls:
        M = int(math.ceil(F*L/math.pi))+8; r = analyse(name, mo, L, M, F, 1e-8, W=wmatrix(mo, L, M), **kw)
        out['fakes'].append(dict(model=name, **{k: r.get(k) for k in ['L', 'M', 'dimS', 'lam_min_W', 'lam_min_perp', 'n_neg_W', 'n_neg_perp']})); print(out['fakes'][-1], flush=True)
fit = {}
for name in ['EPS', 'FTH', 'DH']:
    rr = [r for r in out['fakes'] if r['model'] == name and r['lam_min_perp'] < 0]
    L = np.array([r['L'] for r in rr]); y = np.log(-np.array([r['lam_min_perp'] for r in rr]))
    b, a = np.polyfit(L, y, 1); fit[name] = dict(form='|lam_min| = exp(a + b L)', a=float(a), b=float(b), max_resid=float(np.abs(a+b*L-y).max()))
for F in [45, 90]:
    rr = [r for r in out['zeta'] if r['F'] == F and r['lam_min_W'] > 1e-12*r['normW']]
    L = np.array([r['L'] for r in rr]); y = np.log(np.array([r['lam_min_W'] for r in rr]))
    b, a = np.polyfit(L, y, 1); fit[f'ZETA_F{F}'] = dict(form='lam_min = exp(a + b L), L above floor only', Ls=L.tolist(), a=float(a), b=float(b), max_resid=float(np.abs(a+b*L-y).max()))
out['fits'] = fit; print(json.dumps(fit, indent=1))
json.dump(out, open('fits.json', 'w'), indent=1)
