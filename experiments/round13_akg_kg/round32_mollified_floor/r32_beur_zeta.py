"""Round 32: evaluate the discrete Beurling zeta zeta_B(w) near w = 1/2 (continuation, no zeta zeros):
  log zeta_B(w) = sum_{p^k <= X} p^{-kw}/k  - int_20^X x^{-w} dPi_C  +  [E1((w-1) log 20) - sum_q E1((w-rho_q) log 20)]
(the bracket is the closed-form continuation of int_20^inf x^{-w} dPi_C, dPi_C = (1 - sum_q x^{rho_q-1}) dx / log x).
The truncation error is O(|w| X^{-Re w}) (|Pi_B - Pi_C| <= 1 CHECKED).  Used to explain the Beurling Moebius-floor numerics."""
import json, math, numpy as np, mpmath
X = 1e8; P = np.fromfile('beur_primes.bin'); q = [0.6+28j, 0.6-28j, 0.4+28j, 0.4-28j]
lp = np.log(P)
u = np.linspace(math.log(20), math.log(X), 400001); du = u[1] - u[0]
xs = np.exp(u); dens = 1 - sum((xs**complex(r - 1)).real for r in q)
def logzB(w):
    s = 0j; k = 1
    while True:
        sel = k*lp <= math.log(X)
        if not sel.any(): break
        s += np.sum(np.exp(-k*w*lp[sel]))/k; k += 1
    f = np.exp((1 - w)*u)*dens/u                       # x^{-w} dPi_C/du = x^{1-w} dens / log x
    integ = du*(f.sum() - 0.5*(f[0] + f[-1]))
    cont = complex(mpmath.e1((w - 1)*math.log(20))) - sum(complex(mpmath.e1((w - r)*math.log(20))) for r in q)
    return s - integ + cont
out = {}
for w in (2.0, 1.5):
    out[f'check_w{w}'] = complex(np.exp(logzB(w))).real
grid = []
for sig in (0.5, 0.55, 0.6, 0.7):
    for t in np.arange(-1.6, 1.61, 0.1):
        z = np.exp(logzB(complex(sig, t))); grid.append((sig, float(t), abs(z), float(np.angle(z))))
out['grid'] = grid
z0 = np.exp(logzB(0.5)); h = 1e-4
z1 = (np.exp(logzB(0.5 + h)) - np.exp(logzB(0.5 - h)))/(2*h)
out['zetaB_half'] = [z0.real, z0.imag]; out['dzetaB_half'] = [z1.real, z1.imag]
out['slope'] = (1/z0).real; out['intercept'] = (-z1/z0**2).real
json.dump(out, open('r32_beur_zeta.json', 'w'), indent=1)
print({k: v for k, v in out.items() if k != 'grid'})
for g in grid:
    if g[2] < 0.05: print('small', g)
print('min |zeta_B| on grid', min(grid, key=lambda g: g[2]))
