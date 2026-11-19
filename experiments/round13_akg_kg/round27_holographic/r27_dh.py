"""Round 27: localized Galerkin at the controls' first off-line heights (sector-resolved).
Basis w(u) {cos, sin}(T_k u), T_k on a grid around the target height; window X = e^L. ZETA run at the same heights as
the positive control. NUMERICAL, prime side only."""
import numpy as np, json, math, sys
import r27_sectors as R
from r22_models import lam_table
jobs = [('EPS', 15.67, (55, 100, 300)), ('DH', 85.70, (300, 1000, 3000)), ('ZETA', 85.70, (300, 1000)), ('ZETA', 15.67, (100,))]
out = {}
for nm, T0, Xs in jobs:
    Ts = T0 + np.arange(-6, 6.01, 1.0)
    def basis(a, M, parity, Ts=Ts):
        w = lambda u: np.where(np.abs(u) <= a, np.cos(np.pi * u / (2 * a))**2, 0.0)
        trig = np.cos if parity > 0 else np.sin
        return lambda u: w(u)[..., None] * trig(Ts * np.asarray(u)[..., None])
    R.basis = basis
    for X in Xs:
        L = math.log(X); lam = lam_table(nm, X + 2); row = {}
        for par, sn in ((1, 'even'), (-1, 'odd')):
            W = R.gram_prime(nm, L, len(Ts), par, lam=lam, R=float(T0 + 200), nr=int(40 * (T0 + 200)) + 1)
            ev = np.linalg.eigvalsh(W / np.abs(W).max())
            row[sn] = dict(lam_min_rel=float(ev[0]), n_neg=int(np.sum(ev < -1e-7)))
        out[f'{nm}_T{T0}_X{X}'] = row; print(nm, T0, X, json.dumps(row), flush=True)
json.dump(out, open('r27_localized.json', 'w'), indent=1)
