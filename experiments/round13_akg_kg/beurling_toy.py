"""Round 13d(b'): a CONSISTENT continuous Beurling prime system with an off-line quartet (DMV 2006 sec.3-5 style, toy).
dPi(x) = 1_{x>=x0} (1 - sum_q x^{rho_q - 1}) dx / log x, quartet rho = 1/2 +- 0.1 +- 28i, x0 = 20.
Then -zeta_C'/zeta_C(s) = int_{x0}^inf x^{-s}(1 - sum_q x^{rho-1}) dx = x0^{1-s}/(s-1) - sum_q x0^{rho-s}/(s-rho):
zeta_C has simple zeros exactly at the quartet, and the prime measure is NONNEGATIVE (checked on a grid).
3-4-1: E = int x^{-sigma}(3+4cos(t log x)+cos(2t log x)) dPsi >= 0 automatically; checked numerically.
"""
import numpy as np, mpmath as mp, json
mp.mp.dps = 20
x0 = 20
q = [mp.mpf('0.6')+28j, mp.mpf('0.6')-28j, mp.mpf('0.4')+28j, mp.mpf('0.4')-28j]
xs = np.geomspace(x0, 1e8, 400000)
dens = 1 - sum((xs**complex(r-1)).real for r in q)
D = lambda s: mp.mpf(x0)**(1-s)/(s-1) - mp.fsum(mp.mpf(x0)**(r-s)/(s-r) for r in q)
E = lambda sig, t: float(3*D(mp.mpf(sig)).real + 4*D(sig+1j*t).real + D(sig+2j*t).real)
grid = [(s, t) for s in (1.001, 1.01, 1.1, 1.3, 1.6, 2.0) for t in [28.0, 14.0] + list(np.linspace(0.5, 60, 50))]
vals = [(E(s, t), s, t) for s, t in grid]
mn = min(vals)
res = dict(density_min=float(dens.min()), density_negative=bool((dens < 0).any()), E_min=mn[0], E_min_at=(mn[1], mn[2]),
           residue_check=[float(mp.re(-(r - r*1.000000001)*D(r*1.000000001))) for r in q[:1]])
json.dump(res, open('beurling_toy_results.json', 'w'), indent=1); print(res)
