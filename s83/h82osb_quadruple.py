"""Grok S83 H82-OS-b probe (NUMERICAL, float64, RH OPEN). Zero-side toy of a Beurling/DMV-type fake: the zeta Weil kernel W_zeta
(S14 normalisation, s82os/os_hankel.py) plus an off-line zero quadruple 1/2 +- delta +- i*gamma, kernel K(r) = 4 cosh(delta r) cos(gamma r)
(delta = 0 would be an extra on-line double zero, which only adds a PSD term). NOT the actual Diamond-Montgomery-Vorhauer system
(its primes were not built); it models only the zero-side feature DMV-type systems have (zeros off the line, positive prime measure).
Reports, for the even/odd sectors of [-L,L]: full Weil min eig, number of negative eigenvalues of the pole-free form A = Q - pole
(n_-(A+), n_-(A-)), and nodality of the lowest new negative direction. Usage: python h82osb_quadruple.py"""
import sys, json, math, os
import numpy as np
from numpy.polynomial import legendre as leg
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's82os'))
from os_hankel import build
def basisP(L, N, x):
    t = 2 * x / L - 1
    return np.array([math.sqrt((2 * k + 1) / L) * leg.legval(t, [0] * k + [1]) for k in range(N)])
def kforms(L, N, gam, dl, M=200):
    xg, wg = leg.leggauss(M); x = L / 2 * xg + L / 2; w = L / 2 * wg; P = basisP(L, N, x) * w
    K = lambda r: 4 * np.cosh(dl * r) * np.cos(gam * r)
    return P @ K(x[:, None] - x[None, :]) @ P.T, P @ K(x[:, None] + x[None, :]) @ P.T
def signchanges(L, N, v):
    xs = np.linspace(0, L, 801); g = basisP(L, N, xs).T @ v; g = g[abs(g) > 1e-9 * abs(g).max()]
    return int((np.diff(np.sign(g)) != 0).sum())
res = []
for gam, dl in ((14.134725, 0.0), (14.134725, 0.1), (14.134725, 0.3), (5.0, 0.0), (5.0, 0.3), (5.0, 1.0)):
    for L in (0.3, 0.5, 0.7, 0.9, 1.1):
        N = 32
        T, B = build('ZETA', L, N); c = build.comp; ape = c['ap'] + c['am']; apo = c['ap'] - c['am']
        TK, BK = kforms(L, N, gam, dl)
        Qe = T + B + TK + BK; Qo = T - B + TK - BK
        Ae = Qe - np.outer(ape, ape); Ao = Qo + np.outer(apo, apo)
        ee = np.linalg.eigvalsh(Qe); eo = np.linalg.eigvalsh(Qo); ae, Va = np.linalg.eigh(Ae); ao = np.linalg.eigvalsh(Ao)
        r = dict(kind='float', gamma=gam, delta=dl, L=L, N=N, even_min=float(ee[0]), odd_min=float(eo[0]),
                 even_neg=int((ee < -1e-10).sum()), odd_neg=int((eo < -1e-10).sum()),
                 nA_plus=int((ae < -1e-10).sum()), nA_minus=int((ao < -1e-10).sum()),
                 A_plus_low3=[float(z) for z in ae[:3]], A_minus_low2=[float(z) for z in ao[:2]],
                 signchanges_A_plus_ground=signchanges(L, N, Va[:, 0]),
                 signchanges_A_plus_2nd=signchanges(L, N, Va[:, 1]))
        if ee[0] < -1e-10:
            _, Ve = np.linalg.eigh(Qe); r['signchanges_even_neg'] = signchanges(L, N, Ve[:, 0])
        res.append(r); print(json.dumps(r), flush=True)
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'h82osb_results.json'), 'w'), indent=1)
