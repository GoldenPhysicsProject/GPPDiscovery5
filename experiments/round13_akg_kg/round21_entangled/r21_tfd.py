"""Round 21, Task 1-2: the Bost-Connes / primon-gas thermofield double with n paired to 1/n.
|Psi_beta> = Z^{-1/2} sum_n sqrt(c(n)) n^{-beta/2} |1/n>_L |n>_R   (left copy sits at x = 1/n < 1, right copy at x = n > 1;
the boundary is x = 1, the J-fixed point).  Schmidt weights p(n) = c(n) n^{-beta}/Z.  For this to be a state one needs c(n) >= 0.
Measures: Schmidt positivity, entanglement entropy S, and the total correlation between prime sectors
TC = sum_p H(v_p) - H(n)  (= 0 iff the TFD is a tensor product over primes, i.e. an Euler product in entanglement language).
Models: ZETA (c=1), EPS (r_Q(n), x^2+5y^2), FTH0.2 (coeffs of zeta(s-th)zeta(s+th)), DH (periodic, signed)."""
import json, os, sys, math, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r21_models import DH_A
import mpmath as mp
t0 = time.time(); N = 1_000_000
n = np.arange(N + 1, dtype=float)
def coeff_eps():
    r = np.zeros(N + 1); m = np.arange(-int(math.isqrt(N)), int(math.isqrt(N)) + 1)
    for k in range(-int(math.isqrt(N // 5)), int(math.isqrt(N // 5)) + 1):
        v = m * m + 5 * k * k; v = v[(v <= N) & (v > 0)]; np.add.at(r, v, 1)
    return r
def coeff_fth(th=0.2):
    c = np.zeros(N + 1)
    for a in range(1, N + 1):
        b = np.arange(1, N // a + 1); c[a * b] += a**th * b**(-th)
    return c
C = {'ZETA': np.r_[0.0, np.ones(N)], 'EPS': coeff_eps(), 'FTH0.2': coeff_fth(), 'DH': np.array([DH_A[k % 5] for k in range(N + 1)])}
C['DH'][0] = 0.0
PR = [2, 3, 5, 7, 11, 13]
def vp(p):
    v = np.zeros(N + 1, int); q = p
    while q <= N: v[::q] += 1; q *= p
    return v
V = {p: vp(p) for p in PR}
def H(p): p = p[p > 0]; return float(-(p * np.log(p)).sum())
out = {}
for beta in (2.0, 3.0):
    for name, c in C.items():
        w = c.copy(); w[1:] *= n[1:]**(-beta)
        neg = float(-w[w < 0].sum() / np.abs(w).sum())
        r = dict(min_coeff=float(c[1:].min()), negative_weight_fraction=neg)
        if neg == 0:
            p = w / w.sum(); S = H(p); r['S_entanglement'] = S
            Hm = {}
            for q in PR:
                marg = np.bincount(V[q], weights=p); Hm[q] = H(marg)
            # total correlation among the sectors (v_2, v_3, v_5, v_7, v_11, v_13, rest)
            rest = n.copy(); 
            for q in PR: rest = rest / q**V[q]
            _, inv = np.unique(rest[1:], return_inverse=True)
            Hrest = H(np.bincount(inv, weights=p[1:]))
            r['TC_prime_sectors'] = sum(Hm.values()) + Hrest - S
            j = np.bincount(V[2] * 64 + V[3], weights=p); px = np.bincount(V[2], weights=p); py = np.bincount(V[3], weights=p)
            r['I(v2;v3)'] = H(px) + H(py) - H(j)
        out[f'{name}_beta{beta}'] = r
        print(name, beta, r, flush=True)
# exact zeta formulas: S(beta) = log zeta - beta zeta'/zeta = sum_p S_p ; Hagedorn divergence at beta = 1 (the pole)
mp.mp.dps = 20
def S_exact(b): z = mp.zeta(b); return mp.log(z) - b * mp.zeta(b, derivative=1) / z
def S_p(p, b): q = mp.mpf(p)**(-b); return -mp.log(1 - q) - q * mp.log(q) / (1 - q)
primes = [int(x) for x in mp.primerange(2, 200000)] if hasattr(mp, 'primerange') else None
if primes is None:
    s = np.ones(200001, bool); s[:2] = False
    for i in range(2, 448):
        if s[i]: s[i*i::i] = False
    primes = [int(v) for v in np.nonzero(s)[0]]
out['ZETA_S_additivity'] = {str(b): dict(S_exact=mp.nstr(S_exact(b), 12), sum_p_le_2e5=mp.nstr(mp.fsum(S_p(p, b) for p in primes), 12))
                            for b in (2, 3)}
out['ZETA_Hagedorn'] = {str(b): mp.nstr(S_exact(mp.mpf(b)), 10) for b in ('2', '1.5', '1.1', '1.01', '1.001')}
out['renyi_note'] = 'S_alpha(beta) = log(zeta(alpha beta)/zeta(beta)^alpha)/(1-alpha): log-singular exactly where zeta(alpha beta)=0, i.e. only at complex Renyi index alpha = rho/beta; never for real alpha >= 1/beta'
out['meta'] = dict(N=N, seconds=round(time.time() - t0, 1))
print(out['ZETA_S_additivity'], out['ZETA_Hagedorn'])
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r21_tfd.json'), 'w'), indent=1)
