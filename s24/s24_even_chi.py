"""S24d (NUMERICAL, Galerkin): even susceptibility chi_e = -2 w^T B_e^{-1} w (B_e index 1); Q_e = B_e + 2ww^T > 0 iff B_e|w-perp > 0 and chi_e > 1."""
import sys, os, json; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s17"))
from s17_build import *
import mpmath as mp
mp.mp.dps = 60
def tomp(X, I): return mp.matrix([[mp.mpf(X[i, j].mid().str(70, radius=False)) for j in I] for i in I])
M = int(sys.argv[1]); out = []
for hs in sys.argv[2:]:
    h = arb(3).log()/2 if hs == "hstar" else arb(hs)
    B = build(M, h); F = forms(B, h); Ie = list(range(0, M, 2))
    Qe, Pe, Ne = tomp(F["zeta"], Ie), tomp(B["P"], Ie), tomp(B["N"], Ie); Be = Qe - Pe
    w, V = mp.eigsy(Pe/2); k = max(range(len(w)), key=lambda i: w[i]); wv = V[:, k]*mp.sqrt(w[k])
    Lc = mp.cholesky(Ne); Li = mp.inverse(Lc); ev = sorted(mp.eigsy(Li*Be*Li.T)[0])
    chi = -2*(wv.T*mp.lu_solve(Be, wv))[0]
    row = {"h": hs, "M": M, "B_even_two_lowest": [mp.nstr(ev[0], 6), mp.nstr(ev[1], 6)], "chi_e": mp.nstr(chi, 12), "chi_e-1": mp.nstr(chi-1, 6),
           "rank1_resid": mp.nstr(mp.norm(Pe/2 - wv*wv.T), 3)}
    out.append(row); print(json.dumps(row), flush=True)
json.dump(out, open(f"out_evenchi_M{M}.json", "w"), indent=1)
