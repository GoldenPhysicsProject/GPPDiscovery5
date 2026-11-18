"""Round 13i: primon gas / Bost-Connes Fock space. Fock-only data (no zeros, no continuation):
(1) residues of the Gram kernel <s|s'> = F(conj(s)+s') at the Plancherel line Re s = 1/2 (delta weight pi*Res);
(2) the Witten-index field  F_mu(u) = sum_n mu_F(n) n^{-1/2} g(log n - u), mu_F = Dirichlet inverse of the
    'number operator' coefficients: Moebius for zeta (Spector's supersymmetric primon gas, Tr(-1)^F n^-s = 1/zeta(s)),
    Beurling Moebius on squarefree Beurling integers for B'', Dirichlet inverse of a = r_Q/2 for Epstein.
    Explicit formula: F_mu(u) = sum_rho G(rho-1/2) e^{(rho-1/2)u} / F'(rho) + O(e^{-u/2}) -> zeros enter as RESONANCES (poles of 1/F)."""
import numpy as np, math, json, mpmath as mp
import harness as h
from akg_kg import make_g
X = h.NMAX
out = {}
# zeta Moebius sieve
mu = np.ones(X+1, dtype=np.int8); mu[0] = 0
isp = np.ones(X+1, bool); isp[:2] = False
for p in range(2, X+1):
    if isp[p]:
        isp[2*p::p] = False; mu[p::p] *= -1; mu[p*p::p*p] = 0
nz = np.nonzero(mu)[0]; nz = nz[nz >= 1]
ZETA = dict(n=nz.astype(float), b=mu[nz].astype(float))
# Epstein: Dirichlet inverse of a = r_Q/2
E = h.epstein_model(); a = E['r']/2.0; a[0] = 0
b = np.zeros(X+1); acc = np.zeros(X+1); ak = np.nonzero(a)[0]; ak = ak[ak >= 2]
for n in range(1, X+1):
    b[n] = 1.0 if n == 1 else -acc[n]
    if b[n] != 0:
        k = ak[:np.searchsorted(ak, X//n, side='right')]
        if len(k): acc[n*k] += a[k]*b[n]
ib = np.nonzero(b)[0]; ib = ib[ib >= 1]
EPS = dict(n=ib.astype(float), b=b[ib])
out['EPS_inverse_check'] = float(sum(a[d]*b[12//d] for d in range(1, 13) if 12 % d == 0))  # should be 0 (n=12)
# B'': primes = real primes < 55 + Beurling primes >= 55; squarefree products <= X with sign (-1)^k
Bm = h.beurling_model(x0=55.0)
bp = [float(p) for p in range(2, 55) if all(p % q for q in range(2, int(p**0.5)+1))]
bp += sorted(set(float(n) for n, c in zip(Bm['n'], Bm['c']) if abs(math.log(n) - c) < 1e-9))
bp = np.array(sorted(bp)); vals = []; sg = []
stack = [(1.0, 0, 1)]
while stack:
    v, i, s = stack.pop(); vals.append(v); sg.append(s)
    lim = X/v; j = i
    while j < len(bp) and bp[j] <= lim:
        stack.append((v*bp[j], j+1, -s)); j += 1
o = np.argsort(vals); BEUR = dict(n=np.array(vals)[o], b=np.array(sg, float)[o])
out['n_squarefree'] = dict(ZETA=int(len(ZETA['n'])), EPS_support=int(len(EPS['n'])), BEUR=int(len(BEUR['n'])))
L = 2.0; g = make_g(L, 'even'); us = np.arange(2.0, math.log(X) - L/2 - 0.05, 0.01)
for name, mo in (('ZETA', ZETA), ('EPS', EPS), ('BEUR', BEUR)):
    lg = np.log(mo['n']); w = mo['b']/np.sqrt(mo['n']); F = np.empty(len(us))
    for k, u in enumerate(us):
        i0, i1 = np.searchsorted(lg, [u-L/2, u+L/2]); F[k] = np.sum(w[i0:i1]*g(lg[i0:i1]-u))
    out[name] = {f'rms_u{a}': float(np.sqrt(np.mean(F[(us >= a) & (us < a+2)]**2))) for a in (2.0, 5.0, 8.0, 11.0) if a+2 <= us[-1]+0.01}
    # crude growth exponent of rms(F) between windows [4,6] and [11,13]
    r1 = np.sqrt(np.mean(F[(us >= 4) & (us < 6)]**2)); r2 = np.sqrt(np.mean(F[(us >= 11) & (us < 13)]**2))
    out[name]['growth_rate_per_u'] = float(math.log(r2/r1)/7.0)
    print(name, out[name], flush=True)
# (1) residues at s = 1 of the three 'partition functions' (delta weight = pi * Res on Re s = 1/2)
out['residue'] = dict(ZETA=1.0, EPS_a_normalised=float(mp.pi/(2*mp.sqrt(5))))
# B'': Res = lim (s-1) zeta_B(s); use Dirichlet-series + analytic tail at s = 1+eps
nB = Bm['n']; cB = Bm['c']; m55 = h.zeta_model()
def logzB(s):
    atoms_n = np.concatenate([m55['n'][m55['n'] < 55], nB]); atoms_c = np.concatenate([m55['c'][m55['n'] < 55], cB])
    ell = math.log(atoms_n.max())
    return np.sum(atoms_c/np.log(atoms_n)*atoms_n**(-s)) + complex(mp.e1((s-1)*ell)) - sum(complex(mp.e1((s-r)*ell)) for r in Bm['quartet'])
out['residue']['BEUR'] = [float(abs(eps*np.exp(logzB(1+eps)))) for eps in (1e-2, 1e-3, 1e-4)]
print(out['residue'], out['EPS_inverse_check'])
json.dump(out, open('primon.json', 'w'), indent=1)
