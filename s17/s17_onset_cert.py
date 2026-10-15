"""CERTIFIED (Arb Rayleigh quotient of rationalized Ritz vector) negative witnesses for the pole/arch pieces of Q_zeta
at prime entry (pieces only; the full Q_zeta is NOT claimed negative). M=32."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s17_build import *
import mpmath as mp
mp.mp.dps = 60
M = 32
tomp = lambda X, I: mp.matrix([[mp.mpf(X[i, j].mid().str(70, radius=False)) for j in I] for i in I])
def cert(Q, N, I):
    L = mp.cholesky(tomp(N, I)); Li = mp.inverse(L); ev, V = mp.eigsy(Li*tomp(Q, I)*Li.T)
    k = min(range(len(ev)), key=lambda i: ev[i]); v = Li.T*V[:, k]
    vq = [arb(fmpq(int(mp.nint(v[i]*mp.mpf(10)**40)), 10**40)) for i in range(len(I))]
    num = sum(vq[a]*vq[b]*Q[I[a], I[b]] for a in range(len(I)) for b in range(len(I)))
    den = sum(vq[a]*vq[b]*N[I[a], I[b]] for a in range(len(I)) for b in range(len(I)))
    return (num/den).str(10)
C0, C1 = consts()
for hs, sec, piece in (("0.38", "odd", "Arch0+Pole"), ("0.40", "odd", "Arch0"), ("0.42", "even", "Arch0+Pole"), ("0.45", "even", "F_theta")):
    h = arb(hs); B = build(M, h); N = B["N"]; arch0 = B["K0"]+N*(Tarb(2*h, 0)-C0)
    Q = {"Arch0": arch0, "Arch0+Pole": arch0+B["P"], "F_theta": forms(B, h)["F_theta"]}[piece]
    I = list(range(0, M, 2)) if sec == "even" else list(range(1, M, 2))
    print(f"h={hs} {sec:4s} {piece:11s}: Q(g)/||g||^2 = {cert(Q, N, I)}  (explicit rational g, Arb)")
