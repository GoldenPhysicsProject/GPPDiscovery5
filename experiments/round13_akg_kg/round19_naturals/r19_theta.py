"""Round 19, N1 + N3. d = 1 frame.
N1: additive self-duality test theta_S(1/x) = sqrt(x) theta_S(x), theta_S(x) = 1 + 2 sum_{k} e^{-pi n_k^2 x}, n_k the
'integers' of a prime system S (all-ones on its multiplicative monoid, i.e. Euler product). Equivalently the d=1 state
F_S(x) = sum_k e^{-pi n_k^2 x^2} - 1/(2x) is shadow-invariant: (J F)(x) = x^{-1} F(1/x) = F(x).
Systems: NAT (rational primes -> N), PERT (2 -> 2.01, rest rational), SPLICE (rational p < 55 + Round-13e Beurling
primes >= 55), DMV20 (Round-13e Beurling primes >= 20 only).
N3: heat flow H_t(z) = int Phi(u) e^{t u^2} cos(z u) du leaves the natural class for t != 0 (ratio H_t/H_0 not constant).
NUMERICAL. No zeros used."""
import sys, os, math, json, numpy as np, mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
T = 30000.0
def primes_upto(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return [float(p) for p in np.nonzero(s)[0]]
def monoid(P, X):
    P = sorted(P); out = [1.0]
    def rec(start, val):
        for i in range(start, len(P)):
            v = val*P[i]
            if v > X: break
            out.append(v); rec(i, v)
    sys.setrecursionlimit(100000); rec(0, 1.0); return np.sort(np.array(out))
import harness as h
def beur_primes(x0):
    B = h.beurling_model(x0=x0); n, c = B['n'], B['c']
    return [float(a) for a, b in zip(n, c) if abs(b - math.log(a)) < 1e-9]
R = primes_upto(int(T)+1)
systems = dict(NAT=R, PERT=[2.01]+R[1:], SPLICE=[p for p in R if p < 55]+[q for q in beur_primes(55.0) if q <= T],
               DMV20=[q for q in beur_primes(20.0) if q <= T])
xs = np.logspace(-7, 0, 57)
out = dict(T=T, x_grid=xs.tolist())
for name, P in systems.items():
    n = monoid(P, T)
    def theta(x):
        return 1 + 2*float(np.sum(np.exp(-math.pi*n*n*x)))
    res = [math.sqrt(x)*theta(x) - theta(1/x) for x in xs]
    first_non_nat = None
    if name != 'NAT':
        nat = np.arange(1, len(n)+1, dtype=float); d = np.nonzero(np.abs(n - nat) > 1e-9)[0]
        first_non_nat = float(n[d[0]]) if len(d) else None
    row = dict(n_integers=int(len(n)), density_N_over_T=len(n)/T, first_integer_differing_from_N=first_non_nat,
               max_abs_residual=float(np.max(np.abs(res))),
               residual_at={f'{x:.1e}': r for x, r in zip(xs[::8], res[::8])})
    out[name] = row; print(name, json.dumps(row), flush=True)
# N3 heat flow
mp.mp.dps = 30
def Phi(u):
    u = abs(u)
    return mp.fsum((2*mp.pi**2*k**4*mp.e**(9*u/2) - 3*mp.pi*k**2*mp.e**(5*u/2))*mp.e**(-mp.pi*k**2*mp.e**(2*u)) for k in range(1, 12))
def H(t, z): return 2*mp.quad(lambda u: Phi(u)*mp.e**(t*u*u)*mp.cos(z*u), [0, 1, 2, 4])
heat = {}
for t in (-0.1, 0.1):
    heat[str(t)] = {str(z): float(H(t, z)/H(0, z)) for z in (0, 4, 8, 12)}
out['heat_ratio_Ht_over_H0'] = heat; print('heat', heat)
json.dump(out, open('r19_theta.json', 'w'), indent=1)
