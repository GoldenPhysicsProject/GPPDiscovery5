"""Round 4 (Turbo 6889159916 item 3): OddFace at the single width L = 0.5. NUMERICAL tests of named steps only.
Kernels exactly as in-tree (KreinHPsi.lean / WeilGramSemilocal.lean): archTerm, primeTermS(S_L), bigPsiI2. No zero data.
Odd face form: kform(K_S,(x,-x),(c/2,-c/2)) = 1/2 sum c_i c_j [K_S(x_i-x_j) - K_S(x_i+x_j)], K_S = pole - prime_S + arch.
(1) Scaling test of C1's margin: x = (k-2.5)h, k=0..5, odd block (3 dims), h in {0.005,...,0.08}.
(2) Spread grid in [0, L/2]: x_i = (L/2) i/N, i=1..N, N = 8 and 12.
(3) a = -K''(0) = odd-face coefficient of m_1^2 (Weil form of psi').
(4) Kink table for prime_S in |d|<=L and crossover T_* where w(t) = 2*sum|a_k|."""
import mpmath as mp, json, time
from archgl import arch_multi
mp.mp.dps = 30
t0 = time.time()
L = mp.mpf('0.5'); Smax = mp.e**(2+L)
def triG(y):
    ay = abs(y)
    if ay <= 1: return mp.mpf(2)/3 - y*y + ay**3/2
    if ay <= 2: return (2-ay)**3/6
    return mp.mpf(0)
def triG2(y):
    ay = abs(y)
    if ay <= 1: return -2 + 3*ay
    if ay <= 2: return 2 - ay
    return mp.mpf(0)
def pp(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return p if m == 1 else None
PP = [(n, pp(n)) for n in range(2, 40) if pp(n) and pp(n) < Smax]
def primeS(d):
    return mp.fsum(mp.log(p)/mp.sqrt(n)*(triG(mp.log(n)+d)+triG(-mp.log(n)+d)) for n,p in PP)
bI2 = 16*mp.sinh(mp.mpf('0.25'))**2
pole = lambda d: 2*bI2**2*mp.cosh(d/2)
hs = ['0.005','0.01','0.02','0.04','0.08']
lags = set()
for h in hs:
    for m in range(6): lags.add(mp.nstr(mp.mpf(m)*mp.mpf(h), 10))
for N in (8, 12):
    for m in range(0, 2*N+1): lags.add(mp.nstr(mp.mpf(m)*(L/2)/N, 25))
lags = sorted(lags, key=lambda s: mp.mpf(s))
av = arch_multi([mp.mpf(s) for s in lags], K=1600)
A = {s: v for s, v in zip(lags, av)}
def key(d):
    for s in lags:
        if abs(mp.mpf(s) - abs(d)) < mp.mpf(10)**-12: return s
    raise KeyError(d)
K = lambda d: pole(d) - primeS(d) + A[key(d)]
out = {"label": "NUMERICAL", "L": "0.5", "S_L_prime_powers": [n for n,_ in PP], "bigPsiI2": mp.nstr(bI2, 15)}
# (1) scaling
rows = []
for h in hs:
    hh = mp.mpf(h); x = [(mp.mpf(k)-mp.mpf('2.5'))*hh for k in range(6)]
    # odd block: pairs (x_k, -x_k) with x_k>0: k=3,4,5 -> y = 0.5h,1.5h,2.5h
    y = [x[k] for k in (3,4,5)]
    G = mp.matrix(3,3); G0 = mp.matrix(3,3)
    for i in range(3):
        for j in range(3):
            G[i,j] = (K(y[i]-y[j]) - K(y[i]+y[j]))/2
            G0[i,j] = G[i,j] + 2*bI2**2*mp.sinh(y[i]/2)*mp.sinh(y[j]/2)   # without the pole burden
    # normalise like C1 (vectors e_k - e_{n-1-k}, norm sqrt2): eigen of 2*G w.r.t. |c|^2 with c_full = (c,-c)
    ev = sorted(mp.eigsy(G, eigvals_only=True)); ev0 = sorted(mp.eigsy(G0, eigvals_only=True))
    # C1 normalisation: c_full unit => c has |c|^2 = 1/2, form = 1/2 c G c *... report eig(G)/1 and also /(h^3 log(1/h)), /h^2
    rows.append(dict(h=h, eig_oddblock=[mp.nstr(v,8) for v in ev], eig_without_pole=[mp.nstr(v,8) for v in ev0],
                     min_over_h3log=mp.nstr(ev[0]/(hh**3*mp.log(1/hh)),6), mid_over_h3log=mp.nstr(ev[1]/(hh**3*mp.log(1/hh)),6),
                     max_over_h2=mp.nstr(ev[2]/hh**2,6)))
out["scaling"] = rows
# (2) spread grids
grid = []
for N in (8, 12):
    y = [(L/2)*i/N for i in range(1, N+1)]
    G = mp.matrix(N,N)
    for i in range(N):
        for j in range(N): G[i,j] = (K(y[i]-y[j]) - K(y[i]+y[j]))/2
    ev = sorted(mp.eigsy(G, eigvals_only=True))
    grid.append(dict(N=N, eigs=[mp.nstr(v,6) for v in ev]))
out["spread_grid"] = grid
# (3) a = -K''(0)
def arch2():
    gl = mp.calculus.quadrature.GaussLegendre(mp.mp); nodes = gl.calc_nodes(6, mp.mp.prec)
    tot = mp.mpf(0); h = mp.pi; lp = mp.log(mp.pi); Kp = 1600
    for k in range(Kp):
        c = (2*k+1)*mp.pi
        for (u, w) in nodes:
            t = c + h*u
            tot += (mp.sin(t/2)/(t/2))**4*t*t*(mp.re(mp.digamma(mp.mpc(mp.mpf('0.25'), t/2))) - lp)*w*h
    T = 2*mp.pi*Kp
    # tail: 16 sin^4(t/2)/t^2 * (log(t/2pi) + O(t^-2)); mean of 16 sin^4 = 6; oscillatory parts O(T^-2)
    tail = 6*(mp.log(T/(2*mp.pi)) + 1)/T
    return -(tot + tail)/mp.pi
a2 = arch2()
p2 = mp.fsum(mp.log(p)/mp.sqrt(n)*2*triG2(mp.log(n)) for n,p in PP)
pl2 = 2*bI2**2/4
Kpp = pl2 - p2 + a2
out["second_derivative_at_0"] = dict(arch2=mp.nstr(a2,10), prime2=mp.nstr(p2,10), pole2=mp.nstr(pl2,10), a_minus_Kpp=mp.nstr(-Kpp,10),
                                     a_without_pole=mp.nstr(-(Kpp-pl2),10))
# (4) kinks of prime_S inside |d|<=L and crossover
kinks = []
for n,p in PP:
    ln = mp.log(n)
    for y0, coef in ((1, -mp.mpf(1)/3), (2, mp.mpf(1)/12)):
        for d in (ln - y0, y0 - ln):
            if abs(d) <= L: kinks.append((mp.nstr(d,6), n, mp.nstr(mp.log(p)/mp.sqrt(n)*coef,6)))
sabs = mp.fsum(abs(mp.mpf(k[2])) for k in kinks)
w = lambda t: mp.re(mp.digamma(mp.mpc(mp.mpf('0.25'), t/2))) - mp.log(mp.pi)
Tst = mp.findroot(lambda t: w(t) - 2*sabs, 100)
out["kinks"] = dict(list=kinks, sum_abs=mp.nstr(sabs,8), T_star=mp.nstr(Tst,8), w_zero=mp.nstr(mp.findroot(w, 7),8),
                    modes_below_T_star=mp.nstr((L/2)*Tst/mp.pi,4))
out["seconds"] = round(time.time()-t0,1)
print(json.dumps(out, indent=1)); json.dump(out, open('oddface_L05_results.json','w'), indent=1)
