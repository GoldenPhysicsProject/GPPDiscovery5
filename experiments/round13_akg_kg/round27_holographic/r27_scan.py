"""Round 27 part (3) scan: lam_min of the Weil form on the J-even and J-odd real sectors (Theta-fixed real form),
for ZETA and the controls FTH, EPS, DH, at window widths L (primes/coefficients n <= e^L). Prime side only (no zeros).
Also: frequency content of the most negative eigenvector (where in t the negativity lives)."""
import numpy as np, json, math, sys
from r22_models import lam_table
from r27_sectors import gram_prime, basis, transforms
Ls = [float(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [1, 2, 3, 4, 5, 6]
names = sys.argv[2].split(',') if len(sys.argv) > 2 else ['ZETA', 'FTH', 'EPS', 'DH']
out = {}
Nmax = int(math.exp(max(Ls))) + 2
lams = {nm: lam_table(nm, Nmax) for nm in names}
for L in Ls:
    M = int(3 * L + 6)
    for nm in names:
        row = {}
        for par in (+1, -1):
            W = gram_prime(nm, L, M, par, lam=lams[nm])
            Wn = W / np.abs(W).max()
            ev, V = np.linalg.eigh(Wn)
            v = V[:, 0]; a = L / 2; B = basis(a, M, par)
            r = np.linspace(0, 120, 2401); C, S = transforms(B, a, r); F = (C if par > 0 else S) @ v
            row['even' if par > 0 else 'odd'] = dict(lam_min_rel=float(ev[0]), n_neg=int(np.sum(ev < -1e-9)),
                                                     peak_t=float(r[np.argmax(np.abs(F))]))
        out[f'{nm}_L{L}'] = row; print(nm, L, M, json.dumps(row), flush=True)
json.dump(out, open(f'r27_scan_{"_".join(names)}_{Ls[0]}-{Ls[-1]}.json', 'w'), indent=1)
