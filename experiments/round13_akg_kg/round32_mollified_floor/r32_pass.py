"""Round 32 (argv: N [thetas]; F_0 = zeta^2 via theta=0), one segmented pass to N (default 1e9) for ZETA, F_theta (theta = 0.1, 0.25): Moebius and von Mangoldt
accumulators at every threshold X (fine L-grid, step H, plus the F13S integer samples X = k^4).
No zeros are inputs. Everything here is finite integer arithmetic in float64.

Per model M with Dirichlet inverse mu_M (mu, or mu_theta(p) = -(p^th + p^-th), mu_theta(p^2) = 1, 0 for p^3 | n):
  A_j(X) = sum_{d<=X} mu_M(d) d^{-1/2} (log d)^j  (j=0,1,2)   -> Moebius floor m1(L) = L A0 - A1, m2 = (L^2 A0 - 2L A1 + A2)/2
  B_j(X) = sum_{d<=X} mu(d) d^{-1} (log d)^j  (j=0,1; ZETA only) -> m^(1)(L) = L B0 - B1
  P_j(X) = sum_{n<=X} Lam_M(n) n^{-1/2} (log n)^j (j=0,1)      -> S(L) = L P0 - P1  (the R31 floor's prime sum)
Also max |mu_M(n)| on dyadic blocks [2^j, 2^{j+1})."""
import sys, math, json, time, numpy as np

N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**9
SEG = 5_000_000
H = 0.005
THS = tuple(float(t) for t in sys.argv[2].split(',')) if len(sys.argv) > 2 else (0.1, 0.25)

Lg = np.arange(H, math.log(N) + 1e-12, H)
kmax = int(round(N**0.25));  kmax -= (kmax**4 > N)
Xk = np.arange(2, kmax + 1, dtype=np.int64)**4
T = np.unique(np.concatenate([np.floor(np.exp(Lg) * (1 + 1e-13)).astype(np.int64), Xk]))   # integer thresholds
T = T[(T >= 1) & (T <= N)]
nT = len(T)

r = int(math.isqrt(N))
sv = np.ones(r + 1, bool); sv[:2] = False
for p in range(2, int(r**0.5) + 1):
    if sv[p]: sv[p*p::p] = False
small = np.nonzero(sv)[0].astype(np.int64)
# prime powers p^k (k >= 1) with p <= sqrt N  (their Lambda is set explicitly; large primes detected by the sieve)
pp_n, pp_l = [], []
for p in small:
    q = int(p)
    while q <= N: pp_n.append(q); pp_l.append(math.log(p)); q *= int(p)
o = np.argsort(pp_n); pp_n = np.array(pp_n, np.int64)[o]; pp_l = np.array(pp_l)[o]

names = ['mu_A0', 'mu_A1', 'mu_A2', 'mu_B0', 'mu_B1', 'lam_P0', 'lam_P1']
for th in THS:
    names += [f'f{th}_A0', f'f{th}_A1', f'f{th}_A2', f'lf{th}_P0', f'lf{th}_P1']
acc = {k: np.zeros(nT + 1) for k in names}
dmax = {k: np.zeros(64) for k in ['mu'] + [f'f{th}' for th in THS]}
fp = {th: (lambda th: lambda p: -(p**th + p**(-th)))(th) for th in THS}

t0 = time.time()
for lo in range(1, N + 1, SEG):
    hi = min(lo + SEG, N + 1)
    n = np.arange(lo, hi, dtype=np.int64)
    prod = np.ones(hi - lo, np.int64); mu = np.ones(hi - lo, np.int8)
    F = {th: np.ones(hi - lo) for th in THS}
    for p in small:
        p = int(p)
        if p >= hi: break
        st = (-lo) % p
        s = slice(st, None, p)
        mu[s] *= -1; prod[s] *= p
        for th in THS: F[th][s] *= fp[th](p)
        p2 = p*p
        if p2 < hi:
            s2 = slice((-lo) % p2, None, p2)
            mu[s2] = 0; prod[s2] *= p
            for th in THS: F[th][s2] /= fp[th](p)
            p3 = p2*p
            if p3 < hi:
                s3 = slice((-lo) % p3, None, p3)
                for th in THS: F[th][s3] = 0.0
    q = n // prod
    big = q > 1
    mu[big] *= -1
    qf = q[big].astype(float)
    for th in THS: F[th][big] *= -(qf**th + qf**(-th))
    nf = n.astype(float); w = nf**-0.5; ln = np.log(nf)
    lam = np.zeros(hi - lo)
    isprime_big = (prod == 1) & (n > r)
    lam[isprime_big] = ln[isprime_big]
    a, b = np.searchsorted(pp_n, lo), np.searchsorted(pp_n, hi)
    lam[pp_n[a:b] - lo] = pp_l[a:b]
    idx = np.searchsorted(T, n, side='left')          # n contributes to thresholds T[j] >= n  -> bin j
    def add(key, wt): acc[key] += np.bincount(idx, weights=wt, minlength=nT + 1)
    muf = mu.astype(float); mw = muf*w
    add('mu_A0', mw); add('mu_A1', mw*ln); add('mu_A2', mw*ln*ln)
    add('mu_B0', muf/nf); add('mu_B1', muf*ln/nf)
    lw = lam*w; add('lam_P0', lw); add('lam_P1', lw*ln)
    for th in THS:
        fw = F[th]*w
        add(f'f{th}_A0', fw); add(f'f{th}_A1', fw*ln); add(f'f{th}_A2', fw*ln*ln)
        lt = lam*(nf**th + nf**-th)*w
        add(f'lf{th}_P0', lt); add(f'lf{th}_P1', lt*ln)
    # dyadic maxima
    j0, j1 = int(math.log2(lo)), int(math.log2(hi - 1))
    for j in range(j0, j1 + 1):
        u, v = max(lo, 2**j) - lo, min(hi, 2**(j + 1)) - lo
        if v <= u: continue
        dmax['mu'][j] = max(dmax['mu'][j], float(np.abs(mu[u:v]).max()))
        for th in THS: dmax[f'f{th}'][j] = max(dmax[f'f{th}'][j], float(np.abs(F[th][u:v]).max()))
    if (lo // SEG) % 20 == 0:
        print(f'{hi-1:.3e}  {time.time()-t0:.0f}s', flush=True)

cum = {k: np.cumsum(v)[:nT].tolist() for k, v in acc.items()}       # value at threshold T[j]
json.dump(dict(N=N, H=H, T=T.tolist(), kmax=kmax, cum=cum, dyadic_max={k: v.tolist() for k, v in dmax.items()}),
          open(f'r32_pass_{N:.0e}'.replace('+', '') + ('' if len(sys.argv) <= 2 else '_th' + sys.argv[2].replace(',', '_')) + '.json', 'w'))
print('done', time.time() - t0, flush=True)
