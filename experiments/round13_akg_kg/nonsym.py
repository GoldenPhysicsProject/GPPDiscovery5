"""Round 13d: non-symmetric arithmetic input (Lambda >= 0, Euler product, 3-4-1) as a source of normalizability.
Toys (HYPOTHETICAL, injected): FE+conjugation-symmetric quartets {1/2 +- b +- i g0}. Zeros never used as construction input.
(a) Landau's lemma on the KG norm: N(sigma,U) = int_2^U e^{-sigma u} Phi_g(u)^2 du (nonnegative integrand) -> abscissa = real singularity.
(b) prime-side positivity of the injected field: fake prime density between integers, and the 3-4-1 combination
    E(sigma,t) = 3D(sigma) + 4 Re D(sigma+it) + Re D(sigma+2it) >= 0, D = -zeta'/zeta (true) or D - sum_q 1/(s-rho_q) (injected).
(c) signed sum: lambda_max(prime - pole) vs lambda_max(prime) (absolute), and lambda_min(arch), from the Galerkin matrices.
"""
import numpy as np, math, json, sys
import mpmath as mp
sys.path.insert(0, '.')
from akg_kg import vonmangoldt_table, make_g, W_parts, Phi, gv
out = {}
# (a)
NS, LAM = vonmangoldt_table(int(math.exp(13.3))+1)
L = 0.5; g = make_g(L, 'even'); W = W_parts(g, L, NS, LAM)['W']
du = 0.004; us = np.arange(2.0, 13.0+du/2, du); F, _ = Phi(g, L, us, NS, LAM)
gt = gv(g, L, np.array([28-0.1j, 28+0.1j]))
T = F + 2*np.real(gt[0]*np.exp(28j*us)*np.exp(0.1*us)) + 2*np.real(gt[1]*np.exp(28j*us)*np.exp(-0.1*us))
rows = []
for sig in (0.05, 0.1, 0.15, 0.2, 0.25, 0.3):
    r = dict(sigma=sig)
    for U in (9.0, 11.0, 13.0):
        m = us <= U; e = np.exp(-sig*us[m])
        r[f'true_U{U:g}'] = float(np.trapezoid(e*F[m]**2, us[m])/W)
        r[f'toy_U{U:g}'] = float(np.trapezoid(e*T[m]**2, us[m])/W)
    rows.append(r)
out['a_landau'] = rows
print('a done', flush=True)
# (b)
mp.mp.dps = 20
D = lambda s: -mp.zeta(s, derivative=1)/mp.zeta(s)
def quartet(b, g0): return [mp.mpf(0.5)+b+1j*g0, mp.mpf(0.5)+b-1j*g0, mp.mpf(0.5)-b+1j*g0, mp.mpf(0.5)-b-1j*g0]
def E(sig, t, q):
    f = lambda s: D(s) - (mp.fsum(1/(s-r) for r in q) if q else 0)
    return float(3*f(mp.mpf(sig)).real + 4*f(sig+1j*t).real + f(sig+2j*t).real)
res = {}
for name, q in (('true', None), ('quartet_b0.1_g28', quartet(mp.mpf('0.1'), 28)), ('quartet_b0.3_g28', quartet(mp.mpf('0.3'), 28)),
                ('quartet_b0.47_g28', quartet(mp.mpf('0.47'), 28)), ('quartet_b0.49_g28', quartet(mp.mpf('0.49'), 28))):
    best = (1e9, None)
    for sig in (1.001, 1.003, 1.01, 1.03, 1.1, 1.3, 1.6):
        for t in [28.0, 27.5, 28.5] + list(np.linspace(1, 60, 40)):
            v = E(sig, t, q)
            if v < best[0]: best = (v, (sig, float(t)))
    res[name] = dict(min_E=best[0], at=best[1]); print(name, res[name], flush=True)
out['b_341'] = res
xs = np.geomspace(10, 1e4, 200000)
dens = {}
for b in (0.1, 0.3):
    d = -(2*xs**(b-0.5)*np.cos(28*np.log(xs)) + 2*xs**(-b-0.5)*np.cos(28*np.log(xs)))
    dens[str(b)] = dict(min=float(d.min()), max=float(d.max()), frac_negative=float(np.mean(d < 0)))
out['b_density'] = dens
# (c)
from defect_kappa import galerkin
cr = []
for L in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0):
    M = int(max(24, 10*L)); d = galerkin(L, M)
    ea = np.linalg.eigvalsh(d['arch'])
    cr.append(dict(L=L, M=M, lam_min_arch=float(ea[0]), lam_max_arch=float(ea[-1]),
                   lam_max_prime=float(np.linalg.eigvalsh(d['prime'])[-1]),
                   lam_max_prime_minus_pole=float(np.linalg.eigvalsh(d['prime']-d['pole'])[-1]),
                   signed_bound=float(ea[0]-np.linalg.eigvalsh(d['prime']-d['pole'])[-1])))
    print(cr[-1], flush=True)
out['c_signed'] = cr
json.dump(out, open('nonsym_results.json', 'w'), indent=1, default=float)
print(json.dumps({k: out[k] for k in ('a_landau', 'b_density')}, indent=1, default=float))
