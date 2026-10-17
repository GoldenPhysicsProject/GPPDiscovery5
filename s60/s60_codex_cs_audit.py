# Grok S60 H60a audit of Codex PR#11 6230675184 items (2),(3).
import mpmath as mp, numpy as np, json
mp.mp.dps = 30
out = {}
def xi(s):
    if s == 0 or s == 1: return mp.mpf(1)/2
    return mp.mpf(1)/2*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
P2 = lambda x: mp.diff(mp.primezeta, x, 2)          # P''(x)=sum_p log^2 p p^{-x}
def BR(s):
    M = int(130/s)+5; v = [P2(m*s) for m in range(1, M+1)]
    return sum(v), sum((m)*v[m] for m in range(1, M))   # R = sum_{m>=2}(m-1)P''(ms)
lx2 = lambda s: mp.diff(lambda z: mp.log(xi(z)), s, 2)
lz2 = lambda s: mp.diff(lambda z: mp.log(mp.zeta(z)), s, 2)
rows = {}
for s in [mp.mpf('1.5'), mp.mpf(2), mp.mpf(3)]:
    B, R = BR(s)
    rhs = B + R - 1/s**2 - 1/(s-1)**2 + mp.psi(1, s/2)/4
    rows[str(s)] = {'B': mp.nstr(B,22), 'R': mp.nstr(R,22), '(logzeta)\'\'': mp.nstr(lz2(s),22),
                    'B+R-(logzeta)\'\'': mp.nstr(B+R-lz2(s),5), '(logxi)\'\'': mp.nstr(lx2(s),22), 'rhs-(logxi)\'\'': mp.nstr(rhs-lx2(s),5)}
out['fisher_identity'] = rows; print(json.dumps(rows, indent=1), flush=True)
# direct prime-sum cross-check of B(3), R(3) (independent of primezeta): primes < 2e7, tail bound int_X^oo log^2x x^{-3}dx*(1.3/log X) crude
N = 20_000_000; sv = np.ones(N, bool); sv[:2] = False
for i in range(2, int(N**0.5)+1):
    if sv[i]: sv[i*i::i] = False
pr = np.nonzero(sv)[0].astype(np.float64)
for s in [3.0]:
    lp = np.log(pr); x = pr**s
    out['direct_prime_B(3)_upto2e7'] = float(np.sum(lp**2/(x-1))); out['direct_prime_R(3)_upto2e7'] = float(np.sum(lp**2/(x-1)**2))
# value at s=1
g, g1 = mp.euler, mp.stieltjes(1)
v_claim = mp.pi**2/8 - 1 - g**2 - 2*g1
out['logxi2_at_1_direct'] = mp.nstr(lx2(mp.mpf(1)), 22)
out['pi2/8-1-gamma^2-2gamma1'] = mp.nstr(v_claim, 22)
out['diff_at_1'] = mp.nstr(lx2(mp.mpf(1)) - v_claim, 5)
out['logxi2_at_0 (symmetry)'] = mp.nstr(lx2(mp.mpf(0)), 22)
out['logxi2_at_1/2'] = mp.nstr(lx2(mp.mpf(1)/2), 22)
# zero-side: (log xi)''(1) = -sum_rho 1/(1-rho)^2 ; check with 2000 zeros pairs (real part)
# (3) A_gauge = diag_p (log p)^4/(p-1)^2
lp = np.log(pr); a = lp**4/(pr-1)**2
m = pr <= 1000
TrA_le, TrA2_le = float(np.sum(a[m])), float(np.sum(a[m]**2))
out['A: TrA(p<=1000)'] = TrA_le; out['A: TrA2(p<=1000)'] = TrA2_le
out['A: TrA2 tail p in (1000,2e7)'] = float(np.sum(a[~m]**2)); out['A: TrA tail p in (1000,2e7)'] = float(np.sum(a[~m]))
X = mp.mpf(N)
tail8 = mp.quad(lambda x: mp.log(x)**8/(x-1)**4, [X-1, mp.inf])   # sum over all n>X-1 bounds prime tail (f decreasing)
tail4 = mp.quad(lambda x: mp.log(x)**4/(x-1)**2, [X-1, mp.inf])
out['A: integer-tail bound n>2e7 of (log)^8/(n-1)^4'] = mp.nstr(tail8,5); out['A: integer-tail bound n>2e7 of (log)^4/(n-1)^2'] = mp.nstr(tail4,5)
TrA = float(np.sum(a)); TrA2 = float(np.sum(a**2))
out['A: TrA (all, tail<)'] = [TrA, float(tail4)]; out['A: TrA2 (all, tail<)'] = [TrA2, float(tail8)]
out['A: ratio TrA2/TrA^2'] = TrA2/TrA**2
out['A: ratio upper using Codex tail 0.411'] = (TrA2_le+0.411)/TrA_le**2
out['A: true tail sum_{p>1000} (logp)^8/(p-1)^4 (<2e7 + bound)'] = float(np.sum(a[~m]**2)) + float(tail8)
out['A: Codex tail bound / true tail'] = 0.411/(float(np.sum(a[~m]**2)) + float(tail8))
# xi side: F(u)=xi(1/2+sqrt u)/xi(1/2)=1+e1 u+e2 u^2+..., e1=sum gamma^-2, p2=e1^2-2e2=sum gamma^-4
x0 = xi(mp.mpf(1)/2); d2 = mp.diff(xi, mp.mpf(1)/2, 2); d4 = mp.diff(xi, mp.mpf(1)/2, 4)
e1 = d2/(2*x0); e2 = d4/(24*x0); p2 = e1**2 - 2*e2
out['xi: e1=sum 1/gamma^2'] = mp.nstr(e1, 15); out['xi: p2=sum 1/gamma^4'] = mp.nstr(p2, 15); out['xi: p2/e1^2'] = mp.nstr(p2/e1**2, 15)
# zero cross-check (first 300 zeros; tail of sum gamma^-4 negligible ~1e-8)
gs = [mp.im(mp.zetazero(k)) for k in range(1, 301)]
out['xi zeros(300): sum gamma^-4'] = mp.nstr(sum(1/gm**4 for gm in gs), 12)
# asymptotic counting: N_A(lam) ~ 2*pi*N_xi(lam) (heuristic, leading order)
for lam in [1e-6, 1e-8]:
    NA = int(np.sum(a > lam)); T = lam**-0.5
    Nxi = T/(2*np.pi)*np.log(T/(2*np.pi*np.e)) + 7/8
    out[f'count lam={lam}'] = {'N_A': NA, 'N_xi(RvM)': Nxi, 'ratio': NA/Nxi}
print(json.dumps(out, indent=1, default=str))
