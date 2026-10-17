"""Grok S58 H58c: lift-wall floor. For f in L^2[-L,L] (real), |fhat|^2 has Fourier support in [-2L,2L], so every invisible
character cos(u t), |u| >= 2L, integrates to 0 against |fhat|^2. Hence for ANY invisible lift m,
   sup_t (comb - m)(t) >= P_max(L) := sup_f  sum_{n<e^{2L}} (2 Lambda(n)/sqrt n) <f, S_{log n} f> / ||f||^2
(shifted versions f e^{i t0 x} plus Kronecker give the same bound for sup over t >= T). This script computes P_max(L)
(and even/odd parts) by a Legendre Ritz basis (lower bounds on P_max, converge from below). NUMERICAL. RH OPEN."""
import numpy as np, sys, json, math, os
PS = [int(x) for x in os.environ.get('PRIMES', '').split(',') if x]  # optional prime-group restriction
from numpy.polynomial.legendre import leggauss, legval
def lam_n(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0
def basis(x, L, M):
    u = x / L; out = np.zeros((M, len(x)))
    for d in range(M):
        c = np.zeros(d + 1); c[d] = 1; out[d] = legval(u, c) * math.sqrt((2 * d + 1) / (2 * L))
    return out
def pmax(L, M=60, NQ=400):
    xg, wg = leggauss(NQ); P = np.zeros((M, M)); comb0 = 0.0
    for n in range(2, int(math.exp(2 * L)) + 2):
        u = math.log(n); w = lam_n(n)
        if w == 0 or u >= 2 * L: continue
        if PS and round(math.exp(w)) not in PS: continue
        c = 2 * w / math.sqrt(n); comb0 += c
        a, b = -L, L - u; x = (b - a) / 2 * xg + (a + b) / 2; ww = (b - a) / 2 * wg
        F = basis(x, L, M); G = basis(x + u, L, M); S = (F * ww) @ G.T
        P += c * (S + S.T) / 2
    ev = np.linalg.eigvalsh(P)
    Pe = P[0::2, 0::2]; Po = P[1::2, 1::2]
    return dict(L=L, M=M, comb0=comb0, Pmax=float(ev[-1]), Pmax_even=float(np.linalg.eigvalsh(Pe)[-1]),
                Pmax_odd=float(np.linalg.eigvalsh(Po)[-1]))
if __name__ == '__main__':
    for Ls in sys.argv[1].split(','):
        for M in (40, 80):
            r = pmax(float(eval(Ls)), M); r['Lstr'] = Ls; r['primes'] = PS or 'all'
            r['Tmin_beta0'] = 2 * math.pi * math.exp(r['Pmax']); print(json.dumps(r), flush=True)
