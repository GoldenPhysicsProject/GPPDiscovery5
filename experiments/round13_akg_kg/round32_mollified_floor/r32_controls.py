"""Round 32: EPS (x^2+5y^2, b_n = r_Q(n)/2) and DH (Davenport-Heilbronn, period-5 coefficients) to N (default 1e7):
their OWN Moebius function mu_M = Dirichlet inverse of b_M, their Lambda_M = (b log) * mu_M, and the same accumulators
as r32_pass.py (A_j for the Moebius floor, P_j for the R31 floor). Also dyadic max |mu_M|. No zeros are inputs."""
import sys, math, json, time, numpy as np
from r22_models import KAP

N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**7
H = 0.005
Lg = np.arange(H, math.log(N) + 1e-12, H)
kmax = int(round(N**0.25)); kmax -= (kmax**4 > N)
Xk = np.arange(2, kmax + 1, dtype=np.int64)**4
T = np.unique(np.concatenate([np.floor(np.exp(Lg)*(1 + 1e-13)).astype(np.int64), Xk])); T = T[(T >= 1) & (T <= N)]

def coeffs(name):
    if name == 'DH':
        return np.array([[0, 1, KAP, -KAP, -1][k % 5] for k in range(N + 1)], float)
    b = np.zeros(N + 1); M = int(math.isqrt(N)) + 1
    for x in range(-M, M + 1):
        y = np.arange(-M, M + 1); v = x*x + 5*y*y; v = v[(v > 0) & (v <= N)]
        np.add.at(b, v, 1)
    return b/2

def inverse(b):
    binv = np.zeros(N + 1); binv[1] = 1.0; acc = np.zeros(N + 1)
    nzb = np.nonzero(b[2:])[0] + 2
    for d in range(1, N//2 + 1):
        if d > 1: binv[d] = -acc[d]
        if binv[d] != 0:
            k = nzb[:np.searchsorted(nzb, N//d, side='right')]
            acc[d*k] += binv[d]*b[k]
    binv[N//2 + 1:] = -acc[N//2 + 1:]
    return binv

out = dict(N=N, H=H, T=T.tolist(), kmax=kmax, cum={}, dyadic_max={})
n = np.arange(N + 1); nf = np.maximum(n, 1).astype(float); w = nf**-0.5; ln = np.log(nf)
idx = np.searchsorted(T, n[1:], side='left')
for name in ('DH', 'EPS'):
    t0 = time.time()
    b = coeffs(name); binv = inverse(b)
    print(name, 'inverse', time.time() - t0, flush=True)
    blog = b*ln; lam = np.zeros(N + 1)
    nzi = np.nonzero(binv)[0]
    for d in np.nonzero(blog)[0]:
        if d < 2: continue
        k = nzi[:np.searchsorted(nzi, N//d, side='right')]
        lam[d*k] += blog[d]*binv[k]
    print(name, 'lambda', time.time() - t0, flush=True)
    def cum(wt): return np.cumsum(np.bincount(idx, weights=wt[1:], minlength=len(T) + 1))[:len(T)].tolist()
    mw = binv*w
    out['cum'][f'{name}_A0'] = cum(mw); out['cum'][f'{name}_A1'] = cum(mw*ln); out['cum'][f'{name}_A2'] = cum(mw*ln*ln)
    out['cum'][f'{name}_P0'] = cum(lam*w); out['cum'][f'{name}_P1'] = cum(lam*w*ln)
    out['dyadic_max'][name] = [float(np.abs(binv[2**j:min(2**(j+1), N+1)]).max()) for j in range(int(math.log2(N)) + 1)]
    out[f'{name}_check'] = dict(b=b[1:13].tolist(), mu=binv[1:13].tolist(), lam=lam[1:13].tolist())
    del b, binv, lam, blog
json.dump(out, open(f'r32_controls_{N:.0e}.json'.replace('+', ''), 'w'))
print('done', flush=True)
