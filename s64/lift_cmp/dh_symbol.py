"""Grok S64 H64b design: size of the DH prime-side symbol at L = 1.70 (cost driver for a one-stroke positive certificate).
S_DH(t) = sum_{n<e^{2L}} 2 b(n) n^{-1/2} cos(t log n), b = Dirichlet coefficients of -(L_DH)'/L_DH (S63 dhritz recurrence, float).
Reports: trivial sum |.|, sampled sup over the torus (lower bound for any lift-free A), Galerkin lam_max of C_DH on L^2[-L,L]
(the floor for ANY lift or Schur-type bound). NUMERICAL."""
import math, json, sys
import numpy as np
from numpy.polynomial.legendre import leggauss
L = float(sys.argv[1]); M = int(sys.argv[2])
s5 = math.sqrt(5); kap = (math.sqrt(10-2*s5)-2)/(s5-1)
chi = {0:(0,0),1:(1,0),2:(0,1),3:(0,-1),4:(-1,0)}
N = int(math.exp(2*L)) + 1
a = {n: chi[n%5][0] + kap*chi[n%5][1] for n in range(1, N+1)}
b = {1: 0.0}
for n in range(2, N+1): b[n] = a[n]*math.log(n) - sum(b[d]*a[n//d] for d in range(2, n) if n % d == 0)
ns = [n for n in range(2, N+1) if math.log(n) < 2*L and abs(b[n]) > 1e-15]
w = {n: b[n]/math.sqrt(n) for n in ns}
rng = np.random.default_rng(1); best = -1e9
ts = rng.uniform(0, 1e7, 2_000_000)
S = sum(2*w[n]*np.cos(ts*math.log(n)) for n in ns); best = float(S.max())
xg, wg = leggauss(M+40)
def V(x):
    y = x/L; out = np.empty((M, len(x))); p0, p1 = np.ones_like(y), y
    for j in range(M):
        if j == 0: pj = p0
        elif j == 1: pj = p1
        else: p0, p1 = p1, ((2*j-1)*y*p1-(j-1)*p0)/j; pj = p1
        out[j] = math.sqrt((2*j+1)/(2*L))*pj
    return out
C = np.zeros((M, M))
for n in ns:
    sh = math.log(n); lo, hi = -L+sh, L
    x = (hi-lo)/2*xg+(hi+lo)/2; ww = (hi-lo)/2*wg
    Sm = (V(x)*ww) @ V(x-sh).T; C += w[n]*(Sm+Sm.T)
ev = np.linalg.eigvalsh(C)
print(json.dumps({"L": L, "M": M, "n_terms": len(ns), "trivial": 2*sum(abs(v) for v in w.values()), "sampled_sup_S": best,
  "lam_max_C_DH_galerkin": float(ev[-1]), "T_if_A=lam_max+0.25": 2*math.pi*math.exp(ev[-1]+0.25),
  "T_if_A=sampled_sup+0.25": 2*math.pi*math.exp(best+0.25)}))
