"""S24c (NUMERICAL, Galerkin): even sector. Q_even = B_even + P_even, P_even = 2 w w^T (w_j = <b_j, cosh(x/2)>).
Reports lambda_min(Q_even), lambda_min(B_even | w-perp), lambda_min(B_even | mean-zero)."""
import sys, os, json; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s17"))
from s17_build import *
import mpmath as mp, numpy as np
from numpy.polynomial import legendre as L
mp.mp.dps = 60
def tomp(X, I): return mp.matrix([[mp.mpf(X[i, j].mid().str(70, radius=False)) for j in I] for i in I])
def lam(X, N): Lc = mp.cholesky(N); Li = mp.inverse(Lc); return min(mp.eigsy(Li*X*Li.T)[0])
def perp(a):
    a = np.array(a, float); a /= np.linalg.norm(a); n = len(a)
    Qm, _ = np.linalg.qr(np.column_stack([a] + [np.eye(n)[:, k] for k in range(n)])); return mp.matrix(Qm[:, 1:n].tolist())
M = int(sys.argv[1]); out = []
for hs in sys.argv[2:]:
    h = arb(3).log()/2 if hs == "hstar" else arb(hs)
    B = build(M, h); F = forms(B, h); Ie = list(range(0, M, 2))
    Qe, Pe, Ne = tomp(F["zeta"], Ie), tomp(B["P"], Ie), tomp(B["N"], Ie); Be = Qe - Pe
    w, V = mp.eigsy(Pe); k = max(range(len(w)), key=lambda i: w[i]); wv = [float(V[i, k]) for i in range(len(Ie))]
    x, ww = L.leggauss(80); mv = [float(np.sum(ww*L.legval(x, [0]*j+[1])*(1-x**2)**4)) for j in Ie]
    Zw, Zm = perp(wv), perp(mv)
    cosang = abs(np.dot(wv, mv))/np.linalg.norm(wv)/np.linalg.norm(mv)
    row = {"h": hs, "M": M, "lam_Q_even": mp.nstr(lam(Qe, Ne), 8), "lam_B_even_wperp": mp.nstr(lam(Zw.T*Be*Zw, Zw.T*Ne*Zw), 8),
           "lam_B_even_mean0": mp.nstr(lam(Zm.T*Be*Zm, Zm.T*Ne*Zm), 8), "coeff_cos_angle_w_mean": cosang}
    out.append(row); print(json.dumps(row), flush=True)
json.dump(out, open(f"out_evenperp_M{M}.json", "w"), indent=1)
