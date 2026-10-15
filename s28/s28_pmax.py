"""Grok S28: the floor of every invisible lift.
For supp f in [-L,L], g = f * f~ :  (1/2pi) int comb(t) |F(t)|^2 dt = <f, P_L f>,
P_L f(x) = sum_{log n < 2L} Lambda(n) n^{-1/2} [f(x + log n) + f(x - log n)]  (restricted to [-L, L]),
and every invisible m (frequencies |u| >= 2L) has (1/2pi) int m |F|^2 = 0.  Hence for ANY lift,
sup_t (comb - m)(t) >= <f, P_L f>/||f||^2,  i.e.  A_lift >= lambda_max(P_L)   (PROVED, elementary).
By LP duality on the Bohr compactification (= prime torus) the optimal lift constant equals the sup of
int comb dmu over probability measures with Fourier support in (-2L, 2L); the densities |F|^2/||F||^2 are such
measures, and Krein/Fejer-Riesz says (heuristically) they are dense among them => A* = lambda_max(P_L) (HYPOTHESIS: strong duality).
Numerics: P_L in an orthonormal Legendre basis, exact shifted-overlap integrals by Gauss-Legendre.  RH OPEN."""
import numpy as np, math, sys, json
from numpy.polynomial.legendre import legval, leggauss

def mangoldt_terms(L):
    out = []
    for n in range(2, int(math.exp(2 * L)) + 2):
        if math.log(n) >= 2 * L: continue
        for p in range(2, n + 1):
            if n % p == 0: break
        m = n
        while m % p == 0: m //= p
        if m == 1: out.append((n, math.log(p)))
    return out

def phi(i, x, L):  # orthonormal Legendre on [-L, L]
    c = np.zeros(i + 1); c[i] = 1
    return legval(x / L, c) * math.sqrt((2 * i + 1) / (2 * L))

def Pmat(L, N, nq=None):
    nq = nq or (N + 8)
    xg, wg = leggauss(nq)
    P = np.zeros((N, N))
    for n, lp in mangoldt_terms(L):
        u = math.log(n); c = lp / math.sqrt(n)
        # S(u)_ij = int_{-L+u}^{L} phi_i(x) phi_j(x-u) dx  (x in window, x-u in window)
        a, b = -L + u, L
        x = (b - a) / 2 * xg + (b + a) / 2; w = (b - a) / 2 * wg
        Fi = np.array([phi(i, x, L) for i in range(N)]); Fj = np.array([phi(j, x - u, L) for j in range(N)])
        S = (Fi * w) @ Fj.T
        P += c * (S + S.T)
    return P

if __name__ == '__main__':
    res = {}
    for L in [0.5, 0.8, 1.0, 1.19, 1.6, 2.0, 2.5, 3.0]:
        vals = []
        for N in (40, 80, 120):
            P = Pmat(L, N); ev, V = np.linalg.eigh(P); vals.append(ev[-1])
        AL = sum(2 * lp / math.sqrt(n) for n, lp in mangoldt_terms(L))
        res[L] = dict(lam_max=vals[-1], conv=vals, A_L=AL, T1=2 * math.pi * math.exp(vals[-1]))
        print(f'L={L}: lambda_max(P_L) = {vals[-1]:.5f} (N=40,80,120: {vals[0]:.5f},{vals[1]:.5f},{vals[2]:.5f})  raw A_L={AL:.4f}  2e^L={2*math.exp(L):.3f}  4e^L/L={4*math.exp(L)/L:.3f}  T1*={2*math.pi*math.exp(vals[-1]):.1f}')
    json.dump(res, open('pmax.json', 'w'), indent=1)
