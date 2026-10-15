"""S17 odd-floor onset (NUMERICAL Ritz, Arb entries): per g-parity sector, lambda_min of
  zeta pieces: Arch0 := K0+(T0(2h)-c0)N ; Arch0+Pole ; full Q_zeta ; and chi3/chi5quad (chi(2)=-1) full Q.
Shows where the prime current is needed (pole debt on odd g, TwoFacePoleSignature) as h passes log2/2."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s17_build import *
import mpmath as mp
mp.mp.dps = 60
M = int(sys.argv[1]); HS = sys.argv[2:]
tomp = lambda X, I: mp.matrix([[mp.mpf(X[i, j].mid().str(70, radius=False)) for j in I] for i in I])
def lmin(Q, N, I):
    L = mp.cholesky(tomp(N, I)); Li = mp.inverse(L); return min(mp.eigsy(Li*tomp(Q, I)*Li.T, eigvals_only=True))
C0, C1 = consts(); out = []
for hs in HS:
    h = arb(hs); B = build(M, h); F = forms(B, h); N = B["N"]
    arch0 = B["K0"]+N*(Tarb(2*h, 0)-C0)
    for sec, I in (("even", list(range(0, M, 2))), ("odd", list(range(1, M, 2)))):
        r = {"h": hs, "M": M, "sector": sec, "Arch0": mp.nstr(lmin(arch0, N, I), 8), "Arch0+Pole": mp.nstr(lmin(arch0+B["P"], N, I), 8),
             "Q_zeta": mp.nstr(lmin(F["zeta"], N, I), 8), "Q_chi3": mp.nstr(lmin(F["chi3 (mod 3, odd)"], N, I), 8),
             "Q_chi5quad": mp.nstr(lmin(F["chi5 quadratic (mod 5, even)"], N, I), 8), "primes_present": sorted(B["T"].keys())}
        out.append(r); print(json.dumps(r), flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"onset_M{M}.json"), "w"), indent=1)
