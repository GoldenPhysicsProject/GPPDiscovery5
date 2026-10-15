"""S24b (NUMERICAL, Galerkin): Krein-ceiling test for the S22/S23 claim-1 LP.
H24: the best pointwise c in 'm >= c, m - sigma_A in span{cos(u xi), u>=2h} + DC' equals
lambda_min of B = Arch - prime (no pole) on MEAN-ZERO g in L^2(-h,h), BOTH parities (Krein extension / Akhiezer factorization).
Compares lambda_min(B_odd), lambda_min(B_even), lambda_min(B_even | int g = 0) per h."""
import sys, os, json; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s17"))
from s17_build import *
import mpmath as mp
from numpy.polynomial import legendre as L
import numpy as np
mp.mp.dps = 60
def tomp(X, I): return mp.matrix([[mp.mpf(X[i, j].mid().str(70, radius=False)) for j in I] for i in I])
def mean_vec(M):  # int_{-1}^1 P_j(u)(1-u^2)^4 du (h factor irrelevant for nullspace)
    x, w = L.leggauss(60); return [float(np.sum(w*L.legval(x, [0]*j+[1])*(1-x**2)**4)) for j in range(M)]
def lam(X, N):
    Lc = mp.cholesky(N); Li = mp.inverse(Lc); return min(mp.eigsy(Li*X*Li.T)[0])
M = int(sys.argv[1]); out = []
for hs in sys.argv[2:]:
    h = arb(3).log()/2 if hs == "hstar" else (arb(hs.split("/")[0])/arb(hs.split("/")[1]) if "/" in hs else arb(hs))
    B = build(M, h); F = forms(B, h)
    Io, Ie = list(range(1, M, 2)), list(range(0, M, 2))
    Bo, No = tomp(F["zeta"], Io) - tomp(B["P"], Io), tomp(B["N"], Io)
    Be, Ne = tomp(F["zeta"], Ie) - tomp(B["P"], Ie), tomp(B["N"], Ie)
    a = np.array([mean_vec(M)[j] for j in Ie]); a /= np.linalg.norm(a)
    # basis of a-perp (coefficient space)
    Qm, _ = np.linalg.qr(np.column_stack([a] + [np.eye(len(Ie))[:, k] for k in range(len(Ie))]))
    Z = mp.matrix(Qm[:, 1:len(Ie)].tolist())
    row = {"h": hs, "M": M, "lam_B_odd": mp.nstr(lam(Bo, No), 8), "lam_B_even": mp.nstr(lam(Be, Ne), 8),
           "lam_B_even_mean0": mp.nstr(lam(Z.T*Be*Z, Z.T*Ne*Z), 8)}
    out.append(row); print(json.dumps(row), flush=True)
json.dump(out, open(f"out_krein_M{M}.json", "w"), indent=1)
