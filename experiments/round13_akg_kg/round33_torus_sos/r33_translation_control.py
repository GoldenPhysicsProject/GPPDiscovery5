"""Round 33 vs the Codex/Turbo equal-Haar-norm translation control (GPPDiscovery2 a8172e5, audited by Turbo, PR #15).
Translate the torus by tau = (z_p = -1 for all p): coefficient mu(n) -> mu(n)(-1)^{Omega(n)} = |mu(n)|. The centering
-c chi_1 is a constant function (translation-invariant scalar), so it is left alone. Haar data (H_c) are identical for
f_c and its translate; the Kronecker-line Fejer norm (nu_L, -> delta at the identity) is not. N = 1e7, exact step integrals."""
import json, math, numpy as np
N = 10**7; c = -0.6847652360899366
mu = np.ones(N + 1, np.int8); mu[0] = 0; sv = np.ones(N + 1, bool); sv[:2] = False
for p in range(2, N + 1):
    if sv[p]:
        sv[p*p::p] = False; mu[p::p] *= -1; mu[p*p::p*p] = 0
n = np.arange(1, N + 1, dtype=float); lg = np.log(n); l1 = np.log1p(1/n); out = {}
for name, coef in (('mu', mu[1:].astype(float)), ('abs_mu_translate', np.abs(mu[1:]).astype(float))):
    a = coef/np.sqrt(n); a[0] -= c; M = np.cumsum(a); cq = np.cumsum(M*M*l1); cr = np.cumsum(M*l1); cd = np.cumsum(a*a)
    rows = []
    for L in (6.0, 9.0, 12.0, 15.0, math.log(N)):
        i = int(np.floor(math.exp(L)*(1 + 1e-13))) - 1; t = L - lg[i]
        Q = cq[i] - M[i]**2*l1[i] + M[i]**2*t; mc = cr[i] - M[i]*l1[i] + M[i]*t
        G = 2*Q + L*M[i]**2 - 2*M[i]*mc; H = L*cd[i]
        rows.append(dict(L=L, Haar_H=float(H), Q=float(Q), G=float(G), R=float(G - H), base_point_m_minus_cL=float(mc)))
    out[name] = rows
for a, b in zip(out['mu'], out['abs_mu_translate']):
    print('L=%.2f  H: %.6f | %.6f   Q: %.4g | %.4g   G: %.4g | %.4g   R/H: %.3f | %.3f   <f,1>: %.4g | %.4g' %
          (a['L'], a['Haar_H'], b['Haar_H'], a['Q'], b['Q'], a['G'], b['G'], a['R']/a['Haar_H'], b['R']/b['Haar_H'], a['base_point_m_minus_cL'], b['base_point_m_minus_cL']))
json.dump(out, open('r33_translation_control.json', 'w'), indent=1)
