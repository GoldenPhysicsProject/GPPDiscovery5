"""Round 27 part (3): split the zero side of the sector Weil form into the chi = +1 / chi = -1 bookkeeping.
For a Klein orbit of rho = 1/2 + i g (g = t + i d), the J-even sector contribution is 2 Re[c c^T] (d = 0: 2 c c^T, c real)
and for d != 0: 4 Re[c c^T] = 4 (a a^T - b b^T), a = Re C(g), b = Im C(g)  -> the only minus sign in the zero side.
F_theta (zeros rho +- theta) used as illustration, zeros of zeta from zeta_zeros_300.json."""
import numpy as np, json
from r27_sectors import basis, transforms, gram_prime
from r22_models import TH
zz = np.array(json.load(open('zeta_zeros_300.json')))
out = {}
for L in (1.0, 2.0):
    M = int(3 * L + 6)
    for par, nm in ((+1, 'even'), (-1, 'odd')):
        a = L / 2; B = basis(a, M, par)
        g = zz - 1j * TH                              # rho + theta = 1/2 + theta + i gamma  <->  g = gamma - i theta
        C, S = transforms(B, a, g); F = C if par > 0 else S
        A, Bm = F.real, F.imag
        Wplus = 4 * A.T @ A; Wminus = 4 * Bm.T @ Bm   # sum over quartets {rho +- theta, conj}
        W = Wplus - Wminus
        Wp = gram_prime('FTH', L, M, par)
        ev, V = np.linalg.eigh(W); v = V[:, 0]
        out[f'FTH_L{L}_{nm}'] = dict(prime_vs_zero_maxdiff=float(np.abs(W - Wp).max()), lam_min=float(ev[0]),
            plus_part=float(v @ Wplus @ v), minus_part=float(v @ Wminus @ v),
            first_quartet_minus=float(4 * (Bm[0] @ v)**2), first_quartet_plus=float(4 * (A[0] @ v)**2))
        print(L, nm, out[f'FTH_L{L}_{nm}'], flush=True)
json.dump(out, open('r27_quartet.json', 'w'), indent=1)
