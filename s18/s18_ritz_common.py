"""Ritz helpers (S17 conventions): mpmath dps 60 on Arb midpoints; certified Arb Rayleigh quotient of the rationalized vector."""
from s18_common import *
import mpmath as mp
mp.mp.dps = 60
def tomp(X, I): return mp.matrix([[mp.mpf(X[i, j].mid().str(70, radius=False)) for j in I] for i in I])
def geig(Q, N, I):
    q, n = tomp(Q, I), tomp(N, I); L = mp.cholesky(n); Li = mp.inverse(L)
    ev, V = mp.eigsy(Li*q*Li.T); return ev, V, Li
def rq(Q, N, I, v):
    vq = [fmpq(int(mp.nint(v[i]*mp.mpf(10)**40)), 10**40) for i in range(len(I))]
    num = sum(arb(vq[a])*arb(vq[b])*Q[I[a], I[b]] for a in range(len(I)) for b in range(len(I)))
    den = sum(arb(vq[a])*arb(vq[b])*N[I[a], I[b]] for a in range(len(I)) for b in range(len(I)))
    return num/den, vq
def ritz_min(Q, N, I, extra=None):
    ev, V, Li = geig(Q, N, I); k = min(range(len(ev)), key=lambda i: ev[i]); v = Li.T*V[:, k]
    c, vq = rq(Q, N, I, v)
    ex = None
    if extra is not None:
        ex = sum(arb(vq[a])*arb(vq[b])*extra[I[a], I[b]] for a in range(len(I)) for b in range(len(I))) / \
             sum(arb(vq[a])*arb(vq[b])*N[I[a], I[b]] for a in range(len(I)) for b in range(len(I)))
    return ev[k], c, ex
SECT = lambda M: (("even", list(range(0, M, 2))), ("odd", list(range(1, M, 2))))
