# S49 (Grok): where the content lives -- u<0 (critical line, u=-x^2): first Laguerre inequality in u.
# (a) planted off-line quartet into the TRUE xi zero product (100k zeros), detection = min_u (-m'(u)) < 0
# (b) detection threshold delta*(gamma) vs heuristic 1/(gamma sqrt(2 B)), B = sum 1/(g'^2-g^2)^2
# (c) DH actual: -m'_DH(u) on the window around its off-line zero 0.808517+85.699348i (direct, mpmath)
# (d) two-mode family F_c: critical c* where -m'_c(u) first goes negative on u>=0
import json, math, numpy as np, mpmath as mp
G = np.loadtxt('/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt', usecols=1)
G2 = G**2
def bg(u):  # true xi background, vectorized in u (array)
    u = np.atleast_1d(u); out = np.empty_like(u)
    for i in range(0, len(u), 200):
        out[i:i+200] = np.sum(1.0/(u[i:i+200, None] + G2[None, :])**2, axis=1)
    return out
def planted_min(gp, d, ngrid=4001):
    up = d*d - gp*gp; half = 2*d*gp
    u = np.linspace(up - 1.2*half, up + 1.2*half, ngrid)
    w = u - (d + 1j*gp)**2
    tot = bg(u) + 2*np.real(1/w**2)
    i = tot.argmin(); return float(tot[i]), float(u[i])
res = {}
gp = 85.7
rows = []
for d in [0.01, 0.05, 0.308]:
    mn, at = planted_min(gp, d)
    rows.append(dict(gamma=gp, delta=d, min=mn, at_u=at, dip_pred=-1/(2*d*d*gp*gp), background_at_up=float(bg(np.array([d*d - gp*gp]))[0])))
    print(rows[-1])
res['planted_85.7'] = rows
def thresh(gp):
    lo, hi = 1e-4, 50.0
    if planted_min(gp, hi)[0] < 0: return None
    for _ in range(40):
        mid = math.sqrt(lo*hi)
        if planted_min(gp, mid, 2001)[0] < 0: lo = mid
        else: hi = mid
    return hi
th = []
for target in [85.7, 1000, 5000, 20000, 50000]:
    j = np.searchsorted(G, target); gp = 0.5*(G[j-1] + G[j])  # mid-gap
    B = float(bg(np.array([-gp*gp]))[0])
    ds = thresh(gp)
    Delta = 2*math.pi/math.log(gp/(2*math.pi))
    th.append(dict(gamma=gp, gap=float(G[j] - G[j-1]), mean_spacing=Delta, delta_star=ds, heuristic=1/(gp*math.sqrt(2*B)),
                   delta_star_over_spacing=ds/Delta if ds else None))
    print(th[-1])
res['threshold_midgap'] = th
# extrapolate to H=3e12 with lattice heuristic delta* ~ sqrt2*Delta/pi (mid-gap in a uniform lattice)
H = 3e12; DH_ = 2*math.pi/math.log(H/(2*math.pi))
res['lattice_heuristic_at_H'] = dict(spacing=DH_, delta_star=math.sqrt(2)*DH_/math.pi)
print('lattice heuristic at H=3e12: spacing', DH_, ' delta*', math.sqrt(2)*DH_/math.pi)
# (d) two-mode critical c*
mp.mp.dps = 30
def mmp_c(c, U):
    a = [(mp.mpf(4)**k + c)/mp.factorial(2*k) for k in range(120)]
    f0 = mp.polyval(a[::-1], U); f1 = mp.polyval([k*a[k] for k in range(1, 120)][::-1], U)
    f2 = mp.polyval([k*(k-1)*a[k] for k in range(2, 120)][::-1], U)
    return (f1**2 - f0*f2)/f0**2
def mmp_c_closed(c, U):  # closed form via s=sqrt(U); stable for large U
    s = mp.sqrt(U)
    f = lambda v: mp.cosh(2*mp.sqrt(v)) + c*mp.cosh(mp.sqrt(v))
    return -mp.diff(lambda v: mp.log(f(v)), U, 2)
def cmin(c):
    best = None
    for U in np.linspace(1, 150, 300):
        v = mmp_c_closed(c, mp.mpf(U))
        if best is None or v < best[0]: best = (v, U)
    return best
lo, hi = 100.0, 1000.0
for _ in range(30):
    mid = math.sqrt(lo*hi)
    if cmin(mid)[0] < 0: hi = mid
    else: lo = mid
cstar = hi; vb = cmin(cstar*1.05)
print('two-mode c* (first negative -m\' on u>=0) ~', cstar, ' log c* =', math.log(cstar), ' min at c*1.05:', vb)
res['two_mode_cstar'] = dict(cstar=cstar, log_cstar=math.log(cstar), min_at_105pct_cstar=[mp.nstr(vb[0], 8), vb[1]],
                             heuristic='crossover sqrt(u)~log c; Var~1/(16u) vs stiffness~3/(8u^{3/2}) => sqrt u>6 => c>~e^6')
# (c) DH actual on u<0 near its off-line zero
k = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
chl = [0, 1, 1j, -1j, -1]; chb = [0, 1, -1j, 1j, -1]
def L(s, ch): return 5**(-s)*sum(ch[r]*mp.zeta(s, mp.mpf(r)/5) for r in range(1, 5))
def Lam(s): return (5/mp.pi)**((s + 1)/2)*mp.gamma((s + 1)/2)*((1 - 1j*k)/2*L(s, chl) + (1 + 1j*k)/2*L(s, chb))
mp.mp.dps = 25
def fu(v): return mp.re(Lam(mp.mpf(0.5) + 1j*mp.sqrt(-v)))  # v<0
d0, g0 = mp.mpf('0.308517'), mp.mpf('85.699348')
uc = float(d0**2 - g0**2); half = float(2*d0*g0)
grid = np.linspace(uc - 1.3*half, uc + 1.3*half, 241)
vals = []
for U in grid:
    U = mp.mpf(U)
    f0 = fu(U); f1 = mp.diff(fu, U, 1); f2 = mp.diff(fu, U, 2)
    vals.append(float((f1**2 - f0*f2)/f0**2))
vals = np.array(vals)
i = vals.argmin()
print('DH actual -m\'(u) on window [', grid[0], grid[-1], ']: min', vals[i], 'at u', grid[i], ' (window centre', uc, ')')
res['DH_actual_u_negative'] = dict(window=[float(grid[0]), float(grid[-1])], centre=uc, min=float(vals[i]), at=float(grid[i]),
                                   n_negative=int((vals < 0).sum()), n=len(vals))
json.dump(res, open('s49_detect.json', 'w'), indent=1, default=str)
