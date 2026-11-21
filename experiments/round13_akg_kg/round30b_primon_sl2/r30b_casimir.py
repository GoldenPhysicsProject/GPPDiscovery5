"""Round 30-b, part C: unitarity of the principal-series Casimir = a Pontryagin-index count.
Dictionary (PROVED, elementary): a zero rho = 1/2 + i gamma of xi_M is a weight Delta = rho with SL(2) Casimir
C = Delta(1-Delta) = 1/4 + gamma^2.  The rep is unitary principal series iff C >= 1/4 real iff gamma real;
complementary series (Delta real in (0,1)) is excluded for zeta by GPPVerify F9 xiE_real_ne_zero.  An off-line quartet
gives a conjugate pair of non-real Casimirs -> a ghost pair in a Krein space (Pontryagin 1944: a self-adjoint operator
in Pi_kappa has at most kappa eigenvalues in the open upper half-plane).
The Gram/GNS test of the dilation (scaling) flow: D_M(L) = sum_rho (1 - cos gamma L)/gamma^2 = sum_rho (1-cos gamma L)/(C_rho - 1/4)
is computed from the integers (explicit formula, r22_fejer.rho at gamma0 = 0) -- no zeros as inputs.
G_ij = D(L_i) + D(L_j) - D(|L_i - L_j|) is the Gram matrix of the vectors (e^{i gamma L}-1)/gamma in L^2(zeros):
PSD for all grids <=> all C_rho >= 1/4 (Schoenberg) <=> RH_M.  kappa(L) = #negative eigenvalues on [0, L] counts the
resolved ghost pairs (Bombieri 2000: = #off-line conjugate pairs for large truncations, finitely-many case).
Reproduces Round 23 (r23_gram.py) with X = 1e6 and adds kappa(L)."""
import json, math, numpy as np
from r22_fejer import rho
from r22_models import lam_table
X = 10**6; Lmax = math.log(X)
out = {}
for name in ('ZETA', 'FTH', 'EPS', 'DH'):
    lam = lam_table(name, X + 1); out[name] = {}
    for h in (0.25, 0.1, 0.05):
        allL = np.arange(0, int(Lmax/h) + 1)*h
        Dv = np.zeros(len(allL))
        for i, L in enumerate(allL[1:], 1):
            Dv[i] = np.pi*L/2*rho(name, math.exp(L), np.array([0.0]), lam)[0]
        I = np.arange(1, len(allL)); grid = allL[1:]
        G = Dv[I][:, None] + Dv[I][None, :] - Dv[np.abs(I[:, None] - I[None, :])]
        kap = []
        for Lc in (3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, Lmax):
            m = int(np.sum(grid <= Lc + 1e-9)); w = np.linalg.eigvalsh(G[:m, :m])
            tol = 1e-9*max(1.0, np.abs(w).max())
            kap.append(dict(L=round(float(Lc), 3), n=m, kappa=int(np.sum(w < -tol)), min_eig=float(w[0])))
        out[name][str(h)] = dict(kappa_by_L=kap, D_at=dict(zip([str(x) for x in (1, 3, 6, 9, 12)],
                                  [float(Dv[int(round(x/h))]) for x in (1, 3, 6, 9, 12)])))
        print(name, h, [ (k['L'], k['kappa']) for k in kap], flush=True)
        json.dump(out, open('r30b_casimir.json', 'w'), indent=1)
