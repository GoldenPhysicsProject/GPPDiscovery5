"""Round 33 (torus sum-of-squares), one segmented pass to N (default 1e9). No zeros are inputs.
For each model M (ZETA mu; F_theta with theta = 0 (= zeta^2), 0.1, 0.25) with coefficients a_n = mu_M(n) n^{-1/2},
and for the psi-side vector a_n = (Lambda(n) - 1) n^{-1/2} (ZETA), accumulate per threshold bin (same thresholds as R32):
  dg  = sum a_n^2                         (Haar diagonal, counted once)
  q   = sum M(n)^2 log(1+1/n)             (M(n) = sum_{k<=n} a_k: exact int_0^{log X} M(e^u)^2 du, up to the last cell)
  r   = sum M(n)   log(1+1/n)             (exact int_0^{log X} M(e^u) du = the floor m(L), up to the last cell)
  M   = M(n) at the bin's last integer (= R32's A0)."""
import sys, math, json, time, numpy as np
N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**9
SEG = 5_000_000; H = 0.005; THS = (0.0, 0.1, 0.25)
Lg = np.arange(H, math.log(N) + 1e-12, H)
T = np.unique(np.floor(np.exp(Lg)*(1 + 1e-13)).astype(np.int64)); T = T[(T >= 1) & (T <= N)]; nT = len(T)
r = int(math.isqrt(N)); sv = np.ones(r + 1, bool); sv[:2] = False
for p in range(2, int(r**0.5) + 1):
    if sv[p]: sv[p*p::p] = False
small = np.nonzero(sv)[0].astype(np.int64)
pp_n, pp_l = [], []
for p in small:
    q = int(p)
    while q <= N: pp_n.append(q); pp_l.append(math.log(p)); q *= int(p)
o = np.argsort(pp_n); pp_n = np.array(pp_n, np.int64)[o]; pp_l = np.array(pp_l)[o]
keys = ['mu'] + [f'f{t}' for t in THS] + ['lam']
acc = {f'{k}_{s}': np.zeros(nT + 1) for k in keys for s in ('dg', 'q', 'r')}
last = {k: np.zeros(nT + 1) for k in keys}; carry = {k: 0.0 for k in keys}
fp = {th: (lambda th: lambda p: -(p**th + p**(-th)))(th) for th in THS}
t0 = time.time()
for lo in range(1, N + 1, SEG):
    hi = min(lo + SEG, N + 1); n = np.arange(lo, hi, dtype=np.int64)
    prod = np.ones(hi - lo, np.int64); mu = np.ones(hi - lo, np.int8); F = {th: np.ones(hi - lo) for th in THS}
    for p in small:
        p = int(p)
        if p >= hi: break
        s = slice((-lo) % p, None, p); mu[s] *= -1; prod[s] *= p
        for th in THS: F[th][s] *= fp[th](p)
        p2 = p*p
        if p2 < hi:
            s2 = slice((-lo) % p2, None, p2); mu[s2] = 0; prod[s2] *= p
            for th in THS: F[th][s2] /= fp[th](p)
            p3 = p2*p
            if p3 < hi:
                s3 = slice((-lo) % p3, None, p3)
                for th in THS: F[th][s3] = 0.0
    q = n // prod; big = q > 1; mu[big] *= -1; qf = q[big].astype(float)
    for th in THS: F[th][big] *= -(qf**th + qf**(-th))
    nf = n.astype(float); w = nf**-0.5; l1 = np.log1p(1.0/nf)
    lam = np.zeros(hi - lo); ib = (prod == 1) & (n > r); lam[ib] = np.log(nf[ib])
    a, b = np.searchsorted(pp_n, lo), np.searchsorted(pp_n, hi); lam[pp_n[a:b] - lo] = pp_l[a:b]
    vec = {'mu': mu.astype(float)*w, 'lam': (lam - 1.0)*w}
    for th in THS: vec[f'f{th}'] = F[th]*w
    idx = np.searchsorted(T, n, side='left')
    for k, v in vec.items():
        Mrun = carry[k] + np.cumsum(v); carry[k] = float(Mrun[-1])
        acc[f'{k}_dg'] += np.bincount(idx, weights=v*v, minlength=nT + 1)
        acc[f'{k}_q'] += np.bincount(idx, weights=Mrun*Mrun*l1, minlength=nT + 1)
        acc[f'{k}_r'] += np.bincount(idx, weights=Mrun*l1, minlength=nT + 1)
        # value of M at the last integer of each bin present in this segment
        ub, li = np.unique(idx[::-1], return_index=True); last[k][ub] = Mrun[::-1][li]
    if (lo // SEG) % 20 == 0: print(f'{hi-1:.3e}  {time.time()-t0:.0f}s', flush=True)
cum = {k: np.cumsum(v)[:nT].tolist() for k, v in acc.items()}
for k in keys:  # bins are never empty for T[j] >= 1 contiguous integers, so 'last' at T[j] is M(T[j])
    cum[f'{k}_M'] = last[k][:nT].tolist()
json.dump(dict(N=N, H=H, T=T.tolist(), cum=cum), open(f'r33_pass_{N:.0e}.json'.replace('+', ''), 'w'))
print('done', time.time() - t0, flush=True)
