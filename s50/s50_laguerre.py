# S50 (Grok) H50a: higher-order u-Laguerre inequalities for f(u)=xi(1/2+sqrt u)/xi(1/2) (Csordas-Varga F_1, genus 0 in u).
# Normalization (Csordas-Varga 1990): L_n(f)(x) = (1/(2n)!) sum_{k=0}^{2n} (-1)^{k+n} C(2n,k) f^(k)(x) f^(2n-k)(x)
#   <=> |f(x+iy)|^2 = sum_n L_n(x) y^{2n}   (exact generating identity, checked symbolically in part 0).
# Zero form: real u-zeros a_j=-g_j^2 give factors (1 + t/(a_j-x)^2), t=y^2; an off-line pair a=(d+ig)^2, abar gives
#   P_x(t) = (D+b^2)^2 + 2(D-b^2) t + t^2,  D=(Re a - x)^2, b=Im a = 2 d g   (times |a|^-4).
# So L_n(x)/f(x)^2 = [t^n] E_x(t) P_x(t)/P_x(0), E_x = prod_j (1 + w_j t), w_j = 1/(a_j-x)^2.
import json, math, numpy as np, mpmath as mp, sympy as sp, time
T0 = time.time()
ZF = '/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt'
G = np.loadtxt(ZF, usecols=1); G2 = G**2
res = {'zeros_file': ZF, 'n_zeros': int(len(G))}
# ---------- part 0: symbolic checks (PROVED identities) ----------
x, y, t, b, D, m = sp.symbols('x y t b D m', real=True)
f = sp.Function('f')
# generating identity on a test polynomial with complex roots
z = sp.symbols('z')
poly = (z**2 + 3*z + 7)*(z - 2)*(z**2 - z + 5)
lhs = sp.expand(poly.subs(z, x + sp.I*y)*poly.subs(z, x - sp.I*y))
ok = True
for n in range(0, 4):
    Ln = sum((-1)**(k + n)*sp.binomial(2*n, k)*sp.diff(poly, z, k)*sp.diff(poly, z, 2*n - k) for k in range(2*n + 1))/sp.factorial(2*n)
    ok &= sp.simplify(sp.expand(lhs).coeff(y, 2*n) - Ln.subs(z, x)) == 0
# pair factor
al, be = sp.symbols('alpha beta', real=True)
pair = sp.expand(((al - x)**2 + (be - y)**2)*((al - x)**2 + (be + y)**2))
pf = sp.expand(((al - x)**2 + be**2)**2 + 2*y**2*((al - x)**2 - be**2) + y**4)
ok_pair = sp.simplify(pair - pf) == 0
# lattice: mid-gap two-sided lattice spacing c: E(t)=prod_{j>=0}(1+t/(c(j+1/2))^2)^2 = cosh^2(pi sqrt t / c)
tt = sp.symbols('tt', positive=True)
ser = sp.series(sp.cosh(sp.pi*sp.sqrt(tt))**2, tt, 0, 6).removeO()
ok_lat = all(sp.simplify(ser.coeff(tt, k) - sp.Rational(1, 2)*(2*sp.pi)**(2*k)/sp.factorial(2*k)) == 0 for k in range(1, 6))
res['symbolic'] = dict(generating_identity_n0to3=bool(ok), pair_factor=bool(ok_pair), lattice_cosh2_coeffs=bool(ok_lat))
print('symbolic', res['symbolic'])
# lattice bands, B = 4 pi^2 eps^2, eps = delta/spacing (u-spacing c = 2 g Delta, b = 2 d g):
#  n=1: B < 8 ; n>=2 (M=2n): B in M(M-1) -+ sqrt(M(M-1)(4M-6))
def lat_band(n):
    if n == 1: return (0.0, math.sqrt(8)/(2*math.pi))
    M = 2*n; r = math.sqrt(M*(M - 1)*(4*M - 6))
    return (math.sqrt(M*(M - 1) - r)/(2*math.pi), math.sqrt(M*(M - 1) + r)/(2*math.pi))
res['lattice_bands_eps'] = {n: lat_band(n) for n in range(1, 13)}
print('lattice bands', res['lattice_bands_eps'])
# ---------- part 1: actual zeros, planted mid-gap pair ----------
K = 8
def esym(xs, c):
    """elementary symmetric e_0..e_K of what_j = c^2/(a_j-x)^2, for array xs (excluding nothing)."""
    E = np.zeros((len(xs), K + 1)); E[:, 0] = 1.0
    for g2 in G2:
        w = (c*c)/(xs + g2)**2
        E[:, 1:] = E[:, 1:] + w[:, None]*E[:, :-1]
    return E
def bands(target, eps_grid, nx=41):
    j = np.searchsorted(G, target); gp = 0.5*(G[j - 1] + G[j]); gap = G[j] - G[j - 1]
    Dl = 2*math.pi/math.log(gp/(2*math.pi)); c = 2*gp*Dl
    xs, meta = [], []
    for e in eps_grid:
        d = e*Dl; al = d*d - gp*gp; be = 2*d*gp
        for s in np.linspace(-1, 1, nx):
            xs.append(al + be*s); meta.append((e, al, be))
    xs = np.array(xs); E = esym(xs, c)
    out = {}
    for n in range(1, K - 1):
        out[n] = []
    det = np.zeros((len(eps_grid), K - 1), bool); atalpha = np.zeros((len(eps_grid), K - 1))
    for i, e in enumerate(eps_grid):
        sl = slice(i*nx, (i + 1)*nx); _, al, be = meta[i*nx]
        Dh = ((al - xs[sl])/c)**2; bh = (be/c)**2
        p0 = (Dh + bh)**2; p1 = 2*(Dh - bh); p2 = np.ones_like(Dh)
        for n in range(1, K - 1):
            cn = p0*E[sl, n] + p1*E[sl, n - 1] + (p2*E[sl, n - 2] if n >= 2 else 0)
            v = cn/p0
            det[i, n - 1] = v.min() < 0
            atalpha[i, n - 1] = v[nx//2]
    return dict(gamma=float(gp), gap=float(gap), mean_spacing=Dl, gap_over_mean=float(gap/Dl)), det, atalpha
eps_grid = np.round(np.arange(0.01, 3.001, 0.01), 4)
res['actual'] = {}
for target in [85.7, 1000, 50000]:
    info, det, at = bands(target, eps_grid)
    bd = {}
    for n in range(1, K - 1):
        es = eps_grid[det[:, n - 1]]
        # contiguous intervals
        iv = []
        for e in es:
            if iv and abs(e - iv[-1][1] - 0.01) < 1e-9: iv[-1][1] = float(e)
            else: iv.append([float(e), float(e)])
        bd[n] = iv
    need = []
    for i, e in enumerate(eps_grid):
        nn = [n for n in range(1, K - 1) if det[i, n - 1]]
        need.append(nn[0] if nn else None)
    cover = max([eps_grid[i] for i in range(len(eps_grid)) if all(need[k] is not None for k in range(i + 1))], default=None)
    info.update(bands_eps=bd, first_undetected_eps_orders_le6=None if cover is None else float(cover),
                delta_half_in_eps=0.5/info['mean_spacing'],
                order_needed_at={str(e): need[int(round(e*100)) - 1] for e in [0.25, 0.5, 0.75, 1.0, 1.5, 2.0]})
    res['actual'][str(target)] = info
    print(json.dumps(info), flush=True)
print('t', time.time() - T0)
json.dump(res, open('s50_laguerre.json', 'w'), indent=1, default=str)
