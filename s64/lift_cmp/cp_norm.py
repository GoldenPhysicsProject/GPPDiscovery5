"""Grok S64 H64a(ii): compare three upper bounds for the same quantity
   lam_max(C_p | L^2[-L,L]),  C_p = sum_{n<e^{2L}} Lambda(n) n^{-1/2} (S_{log n} + S_{-log n}),  <f,C_p f> = (1/2pi) int S(t)|Ff(t)|^2 dt,
   S(t) = sum 2 Lambda(n) n^{-1/2} cos(t log n)   (the prime term of Q is -<f,C_p f>; same symbol as our lift's S_lift - h).
 (a) Galerkin Legendre lower estimate of lam_max (Rayleigh-Ritz => converges from below; NUMERICAL);
 (b) Mori VIII weighted Schur test ||C_p|| <= sup_x R_w(x), w = 1 + a (x/L)^2, optimised over a on a fine grid (NUMERICAL sup);
 (c) our invisible-harmonic lift A_cert (B&B-certified, from torus_*.json).
usage: cp_norm.py L M"""
import sys, json, math
import numpy as np
from numpy.polynomial.legendre import leggauss, legval
L = float(eval(sys.argv[1])); M = int(sys.argv[2])
def mang(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
ns = [n for n in range(2, int(math.exp(2*L))+1) if mang(n) > 0 and math.log(n) < 2*L]
wts = {n: mang(n)/math.sqrt(n) for n in ns}
# (a) Galerkin with orthonormal Legendre v_j(x) = sqrt((2j+1)/(2L)) P_j(x/L)
xg, wg = leggauss(M + 40)
def V(x):  # (M, len(x))
    y = x / L; out = np.empty((M, len(x)))
    p0, p1 = np.ones_like(y), y
    for j in range(M):
        if j == 0: pj = p0
        elif j == 1: pj = p1
        else: p0, p1 = p1, ((2*j-1)*y*p1 - (j-1)*p0)/j; pj = p1
        out[j] = math.sqrt((2*j+1)/(2*L)) * pj
    return out
C = np.zeros((M, M))
for n in ns:
    a = math.log(n); lo, hi = -L + a, L     # x in I and x - a in I
    if hi <= lo: continue
    x = (hi-lo)/2*xg + (hi+lo)/2; w = (hi-lo)/2*wg
    Vx, Vxa = V(x), V(x - a)
    S = (Vx * w) @ Vxa.T          # <v_i, S_a v_j>
    C += wts[n] * (S + S.T)
ev = np.linalg.eigvalsh(C)
# (b) Schur test with w = 1 + c (x/L)^2 (even weight; R_w even in x): sup over fine grid of x in [0,L]
xs = np.linspace(-L, L, 400001)
def Rw(c):
    w = lambda x: 1 + c*(x/L)**2
    R = np.zeros_like(xs)
    for n in ns:
        a = math.log(n)
        for s in (a, -a):
            y = xs + s; ind = np.abs(y) <= L
            R[ind] += wts[n] * w(y[ind])
    return (R / w(xs)).max()
cs = np.linspace(0, 3, 301); rs = [Rw(c) for c in cs]; ib = int(np.argmin(rs))
out = {"L": sys.argv[1], "M": M, "lam_max_galerkin": float(ev[-1]), "lam_min_galerkin": float(ev[0]),
       "trivial_sum_2Lambda_over_sqrt_n": 2*sum(wts.values()),
       "schur_w1": float(Rw(0.0)), "schur_mori_c2/3": float(Rw(2/3)), "schur_best_c": float(cs[ib]), "schur_best": float(rs[ib])}
print(json.dumps(out))
