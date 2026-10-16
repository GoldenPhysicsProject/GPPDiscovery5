"""S48 H48a(4): Davenport-Heilbronn analogue of G.  f = (1-ik)/2 L(s,chi) + (1+ik)/2 L(s,chibar), q=5, k=(sqrt(10-2sqrt5)-2)/(sqrt5-1)
(phi enters ONLY via the q=5 input: INPUT-CHOICE).  b(n) = coefficients of -f'/f (real).  psi_f(x)=sum_{n<=x} b(n) (no main term: f entire).
G_f(X) = int_1^X psi_f(t) t^{-3/2} dt = 2 sum_{n<=X} b(n)(n^{-1/2}-X^{-1/2}).
Mellin: (1/s) (-f'/f)(s+1/2)/(s+1/2): poles at s=0 (const -2 f'/f(1/2)) and s = rho_f - 1/2 (amplitude X^{beta-1/2}/|rho(rho-1/2)|).
Checks: (i) no real zero of f on [0.5, 4] (NUMERICAL sign scan) -> Landau forbids BOTH one-sided bounds G_f >= -C X^eps, G_f <= C X^eps.
(ii) G_f numerically oscillates with growing two-sided envelope."""
import numpy as np, mpmath as mp, json, math, sys
N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 2_000_000
k = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
chi = np.array([0, 1, 1j, -1j, -1])
n = np.arange(N + 1); a = ((1 - 1j*k)/2*chi[n % 5] + (1 + 1j*k)/2*np.conj(chi[n % 5])).real; a[0] = 0
logn = np.log(np.maximum(n, 1))
# c = a log - sum_{d|n,1<d<n} c(d) a(n/d): accumulate forward
acc = np.zeros(N + 1); c = np.zeros(N + 1)
for d in range(2, N + 1):
    c[d] = a[d]*logn[d] - acc[d]
    if c[d] != 0 and 2*d <= N:
        m = N // d
        acc[2*d::d][:m - 1] += c[d]*a[2:m + 1]
w = np.zeros(N + 1); w[1:] = c[1:]/np.sqrt(n[1:])
cw = np.cumsum(w); cps = np.cumsum(c)
Xs = np.unique(np.round(np.logspace(1, math.log10(N), 60)).astype(int))
Gf = [2*cw[X] - 2*cps[X]/math.sqrt(X) for X in Xs]
mp.mp.dps = 20
def L(s, ch): return 5**(-s)*sum(ch[r]*mp.zeta(s, mp.mpf(r)/5) for r in range(1, 5))
chl = [0, 1, 1j, -1j, -1]; chb = [0, 1, -1j, 1j, -1]
def f(s): return ((1 - 1j*k)/2*L(s, chl) + (1 + 1j*k)/2*L(s, chb))
sig = [0.5 + 0.05*i for i in range(71)]
fv = [float(mp.re(f(s))) for s in sig]
c0f = float(-2*mp.re(mp.diff(f, 0.5)/f(0.5)))
rho = mp.mpf('0.808517') + 85.699348j
amp = float(2/abs(rho*(rho - 0.5)))
out = dict(N=N, X=Xs.tolist(), G_f=Gf, f_real_axis_min=min(fv), f_real_axis_signs=sorted(set(np.sign(fv).tolist())),
           const_pole_s0=c0f, dh_zero=str(rho), amp_coeff=amp,
           note="G_f(X) ~ c0f + amp*X^{0.3085} cos(85.699 log X + phase) + other zeros; zeta envelope 0.0462")
for X, G in zip(Xs[::6], Gf[::6]): print(f"X={X:>9d} G_f={G: .5f}  model amp*X^0.3085={amp*X**0.308517:.5f}")
print('f(sigma) on [0.5,4]: min', min(fv), 'signs', out['f_real_axis_signs'], ' const(s=0)=', c0f)
# crossing points of DH envelope with zeta RH envelope 0.0462 and with sqrt(log X)
for tgt, name in [(0.0462, 'RH envelope'), (None, 'sqrt(log X)')]:
    if tgt: X = (tgt/amp)**(1/0.308517)
    else:
        X = 10.0
        for _ in range(200): X = (math.sqrt(math.log(X))/amp)**(1/0.308517)
    print(f'single DH zero envelope crosses {name} at X ~ {X:.3e}'); out['cross_' + name] = X
json.dump(out, open('s48_dh.json', 'w'), indent=1)
