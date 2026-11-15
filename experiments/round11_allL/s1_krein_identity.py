"""Round 11, script 1. Explicit formula as a Krein (tau-)Gram identity, all L at once.
(a) SYMBOLIC (sympy): kform(pole) = 2 Psi(i/2)^2 (C^2 - S^2), C = sum c cosh(x/2), S = sum c sinh(x/2).
    => on ker l_c the pole pair {Delta=0, Delta=1} (an off-line tau-pair) is exactly the negative square -2Psi^2 S^2.
(b) NUMERICAL test of a specific step: kform(K)(c) = sum_gamma PsiHat(gamma) |sum_i c_i e^{i gamma x_i}|^2
    (the tau-Gram over zeros; under RH a manifest ||T c||^2), K = pole - prime_S + arch exactly as in-tree.
    Zeros are used ONLY as a comparison target (first NZ zeros, tail bounded), never as a proof input.
Kernels: as experiments/round4_turbo/oddface_L05.py (PsiHat(t) = (sin(t/2)/(t/2))^4, Psi = cubic B-spline triG)."""
import numpy as np, mpmath as mp, sympy as sp, json, sys
from scipy import integrate, special
out = {}
# (a) symbolic
n = 3
c = sp.symbols('c0:%d' % n, real=True); x = sp.symbols('x0:%d' % n, real=True); P = sp.symbols('P', positive=True)
kpole = sum(c[i]*c[j]*2*P*sp.cosh((x[i]-x[j])/2) for i in range(n) for j in range(n))
C = sum(c[i]*sp.cosh(x[i]/2) for i in range(n)); S = sum(c[i]*sp.sinh(x[i]/2) for i in range(n))
diff = sp.simplify(sp.expand((kpole - 2*P*(C**2 - S**2)).rewrite(sp.exp)))
out['a_pole_identity_residual'] = str(diff)   # expect 0  (P stands for Psi(i/2)^2 = PsiHat(i/2))
# (b) numerics
PsiHat = lambda t: (np.sinc(t/(2*np.pi)))**4      # (sin(t/2)/(t/2))^4
bI2 = 16*np.sinh(0.25)**2; PH_i2 = bI2**2          # PsiHat(i/2) = 256 sinh^4(1/4)
def triG(y):
    ay = abs(y)
    return 2/3 - y*y + ay**3/2 if ay <= 1 else ((2-ay)**3/6 if ay <= 2 else 0.0)
def pp(m):
    for p in range(2, m+1):
        if m % p == 0:
            k = m
            while k % p == 0: k //= p
            return p if k == 1 else None
def primeS(d, L):
    tot = 0.0
    for m in range(2, int(np.exp(2+L))+2):
        p = pp(m)
        if p: tot += np.log(p)/np.sqrt(m)*(triG(np.log(m)+d) + triG(-np.log(m)+d))
    return tot
w = lambda t: special.digamma(0.25 + 0.5j*t).real - np.log(np.pi)
def arch(d, T=1200.0):
    f = lambda t: PsiHat(t)*w(t)*np.cos(t*d)
    edges = np.arange(0, T+1e-9, 2*np.pi)
    s = sum(integrate.quad(f, a, b, limit=200, epsabs=1e-14)[0] for a, b in zip(edges[:-1], edges[1:]))
    return s/np.pi   # tail beyond T: |.| < 16*log(T)/(pi*3*T^3) ~ 1.2e-9
pole = lambda d: 2*PH_i2*np.cosh(d/2)
NZ = int(sys.argv[1]) if len(sys.argv) > 1 else 600
gam = np.array([float(mp.zetazero(k).imag) for k in range(1, NZ+1)])
GT = gam[-1]
tail_bound = 2*16/(GT**4)*(np.log(GT/(2*np.pi))+1)*GT/(2*np.pi)  # crude: sum_{gamma>GT} PsiHat <= ...
def zero_side(x, c):
    F = (c[None, :]*np.exp(1j*gam[:, None]*x[None, :])).sum(1)
    return float(2*np.sum(PsiHat(gam)*np.abs(F)**2))   # +gamma and -gamma
rng = np.random.default_rng(11)
rows = []
for L in (0.5, 1.0, 2.0, 3.0, 4.5):
    for nn, ker in ((4, True), (6, True), (5, False)):
        xs = np.sort(rng.uniform(0, L, nn)); xs[0], xs[-1] = 0.0, L
        cc = rng.normal(size=nn)
        if ker:
            v = np.cosh(xs/2); cc = cc - v*(cc@v)/(v@v)
        lags = {}
        Kmat = np.zeros((nn, nn))
        for i in range(nn):
            for j in range(nn):
                d = round(abs(xs[i]-xs[j]), 14)
                if d not in lags: lags[d] = pole(d) - primeS(d, L) + arch(d)
                Kmat[i, j] = lags[d]
        lhs = float(cc@Kmat@cc); rhs = zero_side(xs, cc)
        Ssum = float(cc@np.sinh(xs/2))
        rows.append(dict(L=L, n=nn, ker_lc=ker, kform_K=lhs, zero_tauGram=rhs, diff=lhs-rhs,
                         pole_part=float(cc@np.array([[pole(xs[i]-xs[j]) for j in range(nn)] for i in range(nn)])@cc),
                         minus2Psi2S2=-2*PH_i2*Ssum**2))
        print(rows[-1], flush=True)
out['b_rows'] = rows; out['NZ'] = NZ; out['gamma_max'] = GT; out['zero_tail_per_unit_c2_bound'] = tail_bound
out['label'] = '(a) PROVED symbolic; (b) NUMERICAL (zeros = comparison only)'
json.dump(out, open('s1_results.json', 'w'), indent=1)
print(out['a_pole_identity_residual'], tail_bound)
