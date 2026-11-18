"""Validation of Round 14B Galerkin/model data (writes check_models.json)."""
import sys, os, math, json, numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
import weil_galerkin as wg, models14b as md, harness as h
out = {}
# 1. closed-form Galerkin == Round-13 harness Galerkin (zeta, L=2, M=24)
Z = md.zeta(); L, M = 2.0, 24
a, p, q = wg.operators(Z, L, M); a2, p2, q2, _ = h.operators(dict(Z, pole0=True), L, M)
out['zeta_vs_harness_maxdiff'] = dict(arch=float(np.abs(a-a2).max()), pole=float(np.abs(p-p2).max()), prime=float(np.abs(q-q2).max()))
# 2. DH coefficients: -f'/f at s=3 from Lf vs mpmath
D = md.dh(); cf = md.dh_coeff()
f = lambda s, d=0: mp.diff(lambda z: 5**(-z)*sum(cf[r]*mp.zeta(z, mp.mpf(r)/5) for r in range(1, 5)), s, d)
s = mp.mpf(3); lhs = -f(s, 1)/f(s); rhs = float(np.sum(D['c']/D['n']**3))
tail_bound = float(sum(abs(D['c'][-1])*1.0/k**3 for k in range(int(D['n'][-1]), int(D['n'][-1])*50)))
out['DH_logderiv_s3'] = dict(mpmath=float(lhs), series=rhs, rel=float(abs(lhs-rhs)/abs(lhs)))
# 3. F_theta Galerkin == polarization 2 B_zeta(e^{-theta u} g, e^{theta u} g)  (zeros not used)
th = 0.2; F = md.fth(th); L, M, Mb = 2.0, 30, 160
x = np.linspace(-L/2, L/2, 20001); B = wg.basis(L, Mb, x)
g = (1 - (2*x/L)**2)**6*(1 + 0.4*x + 0.3*np.cos(5*x))
coef = lambda fun: np.trapezoid(B*fun[None, :], x, axis=1)
al, ap, am = coef(g), coef(g*np.exp(th*x)), coef(g*np.exp(-th*x))
WZ = sum(wg.operators(Z, L, Mb)[:2]) - wg.operators(Z, L, Mb)[2]
aF, pF, qF = wg.operators(F, L, Mb); WF = aF + pF - qF
out['Ftheta_polarization'] = dict(direct=float(al@WF@al), polar=float(2*am@WZ@ap))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'check_models.json'), 'w'), indent=1)
