"""S22 (NUMERICAL, Ritz/Galerkin): zeta odd sector on the first prime window.
Q_odd = B - 2 v v^T, B = Arch - prime (no pole), v_j = <b_j, sinh(x/2)>.  Galerkin criterion (GPPVerify OddCollisionTrace):
Q_odd >= 0  <=>  B > 0 and chi := 2 v^T B^{-1} v <= 1.  Reports lambda_min(Q_odd), lambda_min(B), chi, per h and M."""
import sys, os, json; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s17"))
from s17_build import *
import mpmath as mp
mp.mp.dps = 60
def tomp(X, I): return mp.matrix([[mp.mpf(X[i, j].mid().str(70, radius=False)) for j in I] for i in I])
M = int(sys.argv[1]); out = []
for hs in sys.argv[2:]:
    h = arb(3).log()/2 if hs == "hstar" else (arb(hs.split("/")[0])/arb(hs.split("/")[1]) if "/" in hs else arb(hs))
    B = build(M, h); F = forms(B, h); I = list(range(1, M, 2))
    Q, P, N = tomp(F["zeta"], I), tomp(B["P"], I), tomp(B["N"], I)
    Bm = Q - P
    L = mp.cholesky(N); Li = mp.inverse(L)
    ev = lambda X: min(mp.eigsy(Li*X*Li.T)[0])
    w, V = mp.eigsy(-P/2); k = max(range(len(w)), key=lambda i: w[i]); v = V[:, k]*mp.sqrt(w[k])
    chi = 2*(v.T*mp.lu_solve(Bm, v))[0]
    row = {"h": hs, "M": M, "lam_Q_odd": mp.nstr(ev(Q), 8), "lam_B_odd": mp.nstr(ev(Bm), 8), "chi": mp.nstr(chi, 12),
           "1-chi": mp.nstr(1-chi, 6), "rank1_resid": mp.nstr(mp.norm(-P/2 - v*v.T), 3)}
    out.append(row); print(json.dumps(row), flush=True)
json.dump(out, open(f"out_susc_M{M}.json", "w"), indent=1)
