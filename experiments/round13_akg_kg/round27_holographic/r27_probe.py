"""Round 27 part (3): sector-resolved frequency probe (R25b-style) on each J-sector of the Theta-fixed real form.
h_T^+(u) = w(u) cos(T u), h_T^-(u) = w(u) sin(T u), w = cos^2(pi u / L) on [-L/2, L/2]. Prime side only.
Reports min_T W(h_T^pm)/||h_T||^2 and where (T) it is attained, for ZETA, FTH, EPS, DH at X = e^L in {100, 300}."""
import numpy as np, json, math, sys
import r27_sectors as R
from r22_models import lam_table
Ts = np.arange(0.0, 100.0001, 0.25)
def mkbasis(par):
    def basis(a, M, parity):
        w = lambda u: np.where(np.abs(u) <= a, np.cos(np.pi * u / (2 * a))**2, 0.0)
        trig = np.cos if parity > 0 else np.sin
        return lambda u: w(u)[..., None] * trig(Ts * np.asarray(u)[..., None])
    return basis
out = {}
names = sys.argv[1].split(',') if len(sys.argv) > 1 else ['ZETA', 'FTH', 'EPS', 'DH']
for X in (100, 300):
    L = math.log(X)
    for nm in names:
        lam = lam_table(nm, int(X) + 2)
        row = {}
        for par in (+1, -1):
            R.basis = mkbasis(par)
            W = R.gram_prime(nm, L, len(Ts), par, lam=lam, R=400.0, nr=16001)
            a = L / 2; x, wq = np.polynomial.legendre.leggauss(2000); u = a * x; wq = a * wq
            Bu = R.basis(a, 0, par)(u); nrm = (Bu**2 * wq[:, None]).sum(0)
            d = np.diag(W) / np.maximum(nrm, 1e-300)
            ok = nrm > 1e-6 * nrm.max()
            i = int(np.argmin(np.where(ok, d, np.inf)))
            row['even' if par > 0 else 'odd'] = dict(min_ratio=float(d[i]), at_T=float(Ts[i]), n_neg_T=int(np.sum((d < 0) & ok)),
                neg_T_ranges=[float(t) for t in Ts[(d < 0) & ok]][:6] + (['...'] if np.sum((d < 0) & ok) > 6 else []))
        out[f'{nm}_X{X}'] = row; print(nm, X, json.dumps(row), flush=True)
json.dump(out, open('r27_probe.json', 'w'), indent=1)
