"""Round 33: EPS and DH to N (default 1e7) with their OWN mu_M (Dirichlet inverse of b_M; construction identical to
R32 r32_controls.py), same accumulators as r33_pass.py (dg, q, r, M) on the same thresholds. No zeros."""
import sys, math, json, time, numpy as np
sys.path.insert(0, '../round32_mollified_floor'); from r22_models import KAP
N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**7
H = 0.005; Lg = np.arange(H, math.log(N) + 1e-12, H)
T = np.unique(np.floor(np.exp(Lg)*(1 + 1e-13)).astype(np.int64)); T = T[(T >= 1) & (T <= N)]
def coeffs(name):
    if name == 'DH': return np.array([[0, 1, KAP, -KAP, -1][k % 5] for k in range(N + 1)], float)
    b = np.zeros(N + 1); M = int(math.isqrt(N)) + 1
    for x in range(-M, M + 1):
        y = np.arange(-M, M + 1); v = x*x + 5*y*y; v = v[(v > 0) & (v <= N)]; np.add.at(b, v, 1)
    return b/2
def inverse(b):
    binv = np.zeros(N + 1); binv[1] = 1.0; acc = np.zeros(N + 1); nzb = np.nonzero(b[2:])[0] + 2
    for d in range(1, N//2 + 1):
        if d > 1: binv[d] = -acc[d]
        if binv[d] != 0:
            k = nzb[:np.searchsorted(nzb, N//d, side='right')]; acc[d*k] += binv[d]*b[k]
    binv[N//2 + 1:] = -acc[N//2 + 1:]
    return binv
n = np.arange(1, N + 1); nf = n.astype(float); w = nf**-0.5; l1 = np.log1p(1/nf)
idx = np.searchsorted(T, n, side='left'); out = dict(N=N, H=H, T=T.tolist(), cum={})
def cum(wt): return np.cumsum(np.bincount(idx, weights=wt, minlength=len(T) + 1))[:len(T)].tolist()
for name in ('DH', 'EPS'):
    t0 = time.time(); a = inverse(coeffs(name))[1:]*w; M = np.cumsum(a)
    out['cum'][f'{name}_dg'] = cum(a*a); out['cum'][f'{name}_q'] = cum(M*M*l1); out['cum'][f'{name}_r'] = cum(M*l1)
    out['cum'][f'{name}_M'] = M[T - 1].tolist(); print(name, time.time() - t0, flush=True)
json.dump(out, open(f'r33_controls_{N:.0e}.json'.replace('+', ''), 'w')); print('done')
