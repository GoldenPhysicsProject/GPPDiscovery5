"""S17 Ritz table: lambda_min(Q_eta)/||g||^2 per g-parity sector, eta in {zeta, F_1/4, chi3, chi4, chi5 quad, chi5 quartic}.
UPPER bounds: Ritz minimum (NUMERICAL, mpmath dps 60 on Arb midpoints) and CERTIFIED (Arb Rayleigh quotient of the
rationalized Ritz vector: an Arb ball for Q(g)/||g||^2 of an explicit g, hence inf <= that).
LOWER bounds (CERTIFIED, Arb, elementary; valid only where the eta prime current vanishes on (-2h,2h)):
  E_eps >= E_comb(h) (S17 bathtub+Schur combined with Astra numerical-radius chains; see s17_bounds.py), P >= 0 on even g, P >= -2(sinh h - h)||g||^2 (Codex/Lean pole bound).
usage: python s17_ritz.py M h1 h2 ...   (h as decimal or a/b)"""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s17_build import *
from s17_bounds import E_comb
import mpmath as mp
mp.mp.dps = 60
M = int(sys.argv[1]); HS = sys.argv[2:]
def harb(s): a, b = (s.split("/")+["1"])[:2]; return arb(a)/arb(b) if "/" in s else arb(s)
def tomp(X, I): return mp.matrix([[mp.mpf(X[i, j].mid().str(70, radius=False)) for j in I] for i in I])
def ritz(Q, N, I):
    q, n = tomp(Q, I), tomp(N, I); L = mp.cholesky(n); Li = mp.inverse(L)
    ev, V = mp.eigsy(Li*q*Li.T); k = min(range(len(ev)), key=lambda i: ev[i]); v = Li.T*V[:, k]
    # certify Rayleigh quotient of rationalized v
    vq = [fmpq(int(mp.nint(v[i]*mp.mpf(10)**40)), 10**40) for i in range(len(I))]
    num = sum(arb(vq[a])*arb(vq[b])*Q[I[a], I[b]] for a in range(len(I)) for b in range(len(I)))
    den = sum(arb(vq[a])*arb(vq[b])*N[I[a], I[b]] for a in range(len(I)) for b in range(len(I)))
    return ev[k], num/den
def lower(nm, h, sector):
    C0, C1 = consts(); pole = 2*(h.sinh()-h)
    if nm == "zeta":
        if arb(2).log() < 2*h: return None
        return E_comb(h, 0)-C0-(0 if sector == "even" else pole)
    if nm == "F_theta":
        if arb(2).log() < 2*h: return None
        th = arb(1)/4
        # int_0^{2h} nu0(t)(cosh th t - 1) dt, entire integrand written as [t nu0(t)]*(th^2 t/2) sinc(i th t/2)^2
        I = acb(0, 1)
        corr = acb.integral(lambda z, an: ((z/2).exp()/(I*z).sinc())*(acb(th)**2*z/2)*((I*acb(th)*z/2).sinc())**2, 0, acb(2*h)).real
        pb = sum(((2*a*h).sinh()/(2*a)-h) for a in (arb(1)/2+th, arb(1)/2-th))
        return 2*(E_comb(h, 0)-C0-corr)-(0 if sector == "even" else 2*pb)
    for cn, f, e, ch in CHARS:
        if cn == nm:
            # prime-free iff Re eta(n)=0 for all prime powers n < e^{2h}
            n = 2
            while arb(n).log() < 2*h:
                if vM(n) and ch(n) != 0: return None
                n += 1
            return E_comb(h, e)+arb(f).log()-(C0 if e == 0 else C1)
rows = []
for hs in HS:
    h = harb(hs); B = build(M, h); F = forms(B, h)
    for nm, Q in F.items():
        for sector, I in (("even", list(range(0, M, 2))), ("odd", list(range(1, M, 2)))):
            lam, cert = ritz(Q, B["N"], I); lb = lower(nm, h, sector)
            row = {"h": hs, "M": M, "eta": nm, "sector": sector, "ritz_min": mp.nstr(lam, 12),
                   "cert_upper": cert.str(12), "cert_lower": (lb.str(12) if lb is not None else "n/a (primes present)")}
            rows.append(row); print(json.dumps(row), flush=True)
json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"ritz_M{M}.json"), "w"), indent=1)
