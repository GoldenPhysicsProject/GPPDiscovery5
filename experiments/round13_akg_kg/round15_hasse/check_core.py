"""B_0 must equal the 14B Weil matrix W = arch + pole - prime (weil_galerkin.operators)."""
import numpy as np, hasse_core as hc, weil_galerkin as wg, json
out = {}
for name in ['ZETA', 'EPS', 'DH', 'FTH']:
    mo = hc.MODELS[name](); l = 0.6; M = hc.Mof(l, 45.0)
    a, p, r = wg.operators(mo, l, M); W = a + p - r
    B = hc.Bs(mo, l, M, 0.0)['B']; B = 0.5*(B + B.T)
    out[name] = dict(maxdiff=float(np.abs(W-B).max()), scale=float(np.abs(W).max()),
                     lam14b=float(np.linalg.eigvalsh(W)[0]), lam15=float(np.linalg.eigvalsh(B)[0]))
    print(name, out[name], flush=True)
# shift consistency: g supported in [-l/2,l/2] and tau_s g inside a bigger window: compare with 14B W on window l+s
mo = hc.MODELS['ZETA'](); l = 0.6; s = 0.9
json.dump(out, open('check_core.json', 'w'), indent=1)
