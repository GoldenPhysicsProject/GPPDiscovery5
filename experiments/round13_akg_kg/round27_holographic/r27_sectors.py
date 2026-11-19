"""Round 27 part (3): the Weil form on the four Klein sectors of ONE plane with two faces.
Faces: J f(x) = x^{-1} f(1/x) (x <-> 1/x, s <-> 1-s). Sides: C f = conj f (s <-> conj s). Theta = CJ (s <-> 1 - conj s).
Log coordinate u = log x, half-width a = L/2: basis b(u) supported in [-a, a] (so g = b*b is supported in [-L, L], n <= e^L).
Sector Gram (proved in ROUND27.md section 3):  W^{+}_{jk} = sum_rho C_j(g_rho) C_k(g_rho),  W^{-}_{jk} = sum_rho S_j S_k,
  C_j(r) = int b_j(u) cos(r u) du,  S_j(r) = int b_j(u) sin(r u) du,  rho = 1/2 + i g_rho. Cross blocks vanish identically.
Prime side: Round 22 explicit formula (r22_models.py), with g_jk = eps * (b_j conv b_k).  No zeros are used on the prime side;
zero side (zeta, F_theta only) is an independent check of the explicit formula. Labels: NUMERICAL unless stated."""
import numpy as np, math, json, sys
from scipy.special import digamma
from r22_models import MODELS, lam_table, TH

def basis(a, M, parity):
    js = np.arange(M) if parity > 0 else np.arange(1, M + 1)
    w = lambda u: np.where(np.abs(u) <= a, np.cos(np.pi * u / (2 * a))**2, 0.0)
    if parity > 0: return lambda u: w(u)[..., None] * np.cos(js * np.pi * np.asarray(u)[..., None] / a)
    return lambda u: w(u)[..., None] * np.sin(js * np.pi * np.asarray(u)[..., None] / a)

def transforms(B, a, r, nq=3000):
    x, wq = np.polynomial.legendre.leggauss(nq); u = a * x; wq = a * wq
    Bu = B(u)                                     # (nq, M)
    r = np.asarray(r)
    if np.iscomplexobj(r):
        Cr = (np.cos(np.multiply.outer(r, u)) * wq) @ Bu; Sr = (np.sin(np.multiply.outer(r, u)) * wq) @ Bu
    else:
        Cr = (np.cos(np.multiply.outer(r, u)) * wq) @ Bu; Sr = (np.sin(np.multiply.outer(r, u)) * wq) @ Bu
    return Cr, Sr

def conv(B, a, xs, nq=1500):
    x, wq = np.polynomial.legendre.leggauss(nq); v = a * x; wq = a * wq
    Bv = B(v)                                     # (nq, M)
    out = np.zeros((len(xs), Bv.shape[1], Bv.shape[1]))
    for i, X in enumerate(xs):
        Bk = B(X - v)                             # (nq, M)
        out[i] = (Bv * wq[:, None]).T @ Bk
    return out

def gram_prime(name, L, M, parity, lam=None, pmax=None, R=600.0, nr=24001):
    a = L / 2; B = basis(a, M, parity); mod = MODELS[name]
    N = int(math.exp(L)) + 1
    if lam is None: lam = lam_table(name, N)
    lam = lam[:N + 1].copy()
    if pmax is not None: lam[int(pmax) + 1:] = 0.0
    ns = np.nonzero(lam)[0]; ns = ns[np.log(ns) < L]
    # h_jk(r) = C_j C_k (even) or S_j S_k (odd); real r grid for arch
    r = np.linspace(-R, R, nr); dr = r[1] - r[0]
    Cr, Sr = transforms(B, a, r); F = Cr if parity > 0 else Sr
    ReGG = sum(c * np.real(digamma(aa + 1j * bb * r)) for (c, aa, bb) in mod['arch'])
    wgt = (mod['logQ'] + ReGG) * dr / np.pi
    W = (F * wgt[:, None]).T @ F
    for z in mod['poles']:
        gz = (z - 0.5) / 1j
        Cz, Sz = transforms(B, a, np.array([gz]))
        Fz = (Cz if parity > 0 else Sz)[0]
        W += 2 * np.real(np.outer(Fz, Fz))
    G = conv(B, a, np.log(ns.astype(float)))     # b_j conv b_k at log n
    eps = 1.0 if parity > 0 else -1.0
    W -= 2 * eps * np.einsum('n,njk->jk', lam[ns] / np.sqrt(ns), G)
    return 0.5 * (W + W.T)

def gram_zero(gammas, L, M, parity):
    a = L / 2; B = basis(a, M, parity)
    g = np.asarray(gammas, complex)
    Cg, Sg = transforms(B, a, g); F = Cg if parity > 0 else Sg
    W = np.real(F.T @ F)                          # sum over listed gammas (list closed under the Klein group)
    return 0.5 * (W + W.T)

def klein(gs):
    out = []
    for gm in gs:   # Klein orbit of rho = 1/2 + i gm: {gm, -gm, conj gm, -conj gm}; collapses to a pair on the line
        orb = {complex(gm), complex(-gm), complex(np.conj(gm)), complex(-np.conj(gm))}
        out += list(orb)
    return out

if __name__ == '__main__':
    zz = json.load(open('zeta_zeros_300.json'))
    res = {}
    # (A) explicit-formula check: prime side vs zero side, zeta and F_theta, both sectors
    for name, gs in (('ZETA', [complex(t) for t in zz]), ('FTH', [complex(t, -TH) for t in zz])):
        for L in (1.0, 2.0):
            for par in (+1, -1):
                Wp = gram_prime(name, L, 8, par); Wz = gram_zero(klein(gs), L, 8, par)
                res[f'check_{name}_L{L}_{"even" if par>0 else "odd"}'] = dict(
                    maxabs_diff=float(np.abs(Wp - Wz).max()), scale=float(np.abs(Wz).max()),
                    lam_min_prime=float(np.linalg.eigvalsh(Wp)[0]), lam_min_zero=float(np.linalg.eigvalsh(Wz)[0]))
                print(name, L, par, res[f'check_{name}_L{L}_{"even" if par>0 else "odd"}'], flush=True)
    json.dump(res, open('r27_sector_check.json', 'w'), indent=1)
