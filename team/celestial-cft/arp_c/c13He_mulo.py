"""C.13.H^e: largest mu_lo with mu_hi*(mu_lo) > c2, and the resulting C.13.H' range (NUMERICAL, Galerkin).
l_H = pi mu'_lo (mu'_hi - c2) / (2 lambda T c2 mu'_hi), mu' = mu - 24.88 l, lambda = 1 + l/a (solved as fixed point).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import numpy as np
from scipy.linalg import eigh
from c13He_even_bandsplit import mats, sector_basis, c2, a
def muhi(Rm, Hlo, Hhi, B, mlo):
    return eigh(B.T@(Rm - mlo*Hlo)@B, B.T@Hhi@B, eigvals_only=True)[0]
def lH(mhi, mlo, T):
    # largest l with l <= f(l) := pi mu'_lo (mu'_hi - c2)/(2 lambda T c2 mu'_hi); f is decreasing in l
    def f(l):
        mh, ml = mhi - 24.88*l, mlo - 24.88*l
        if ml <= 0 or mh <= c2: return -1.0
        return np.pi*ml*(mh - c2)/(2*(1 + l/a)*T*c2*mh)
    lo, hi = 0.0, 0.1
    if f(0) <= 0: return 0.0
    for _ in range(60):
        mid = (lo+hi)/2; (lo, hi) = (mid, hi) if mid <= f(mid) else (lo, mid)
    return lo
for s in ('even', 'odd'):
    for T in (12.875, 14.0, 16.0):
        for M in (320,):
            Rm, H, Hlo, Hhi = mats(M, T); B = sector_basis(M-1, s)
            lo, hi = 0.0, 2.0
            for _ in range(40):
                mid = (lo+hi)/2
                (lo, hi) = (mid, hi) if muhi(Rm, Hlo, Hhi, B, mid) > c2 else (lo, mid)
            best = max(((ml, muhi(Rm, Hlo, Hhi, B, ml)) for ml in np.linspace(0, lo, 81)[1:]), key=lambda p: lH(p[1], p[0], T))
            print(f"{s} T={T} M={M}: max mu_lo with mu_hi*>c2: {lo:.5f}; best pair for C.13.H' (mu_lo, mu_hi)=({best[0]:.5f}, {best[1]:.4f}) -> l_H = {lH(best[1], best[0], T):.2e}", flush=True)
