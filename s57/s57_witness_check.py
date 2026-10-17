"""Grok S57 H57b(i): INDEPENDENT second evaluation of a knockout witness (mpmath only; no flint, no S14 moment code).
Reads the rational coefficients c_j from a s57_witness.py json, builds f(x)=sum c_j P_j(x/L)(1-(x/L)^2)^k (mpmath legendre ->
monomials), computes the autocorrelation g(x)=int f(y)f(y-x)dy on [0,2L] by fixed Gauss-Legendre (exact for the polynomial
integrand), then the Weil functional directly in real space with tanh-sinh quadrature:
  Q = 2 int_0^{2L} g(x) 2cosh(x/2) dx + (-gamma-log pi) g(0) + int_0^{2L} (g0-g)e^{-x}/sinh x dx - g0 log(1-e^{-4L})
      - 2 int_0^{2L} g e^{-x/4} sinh(3x/4)/sinh x dx - sum_{n<e^{2L}, n not a power of p} (2 Lambda(n)/sqrt n) g(log n).
(The archimedean formula is checked against the digamma/Fourier form and the full explicit formula in s57_archcheck.py.)
usage: python s57_witness_check.py witness.json [dps]"""
import sys, json, time
from fractions import Fraction
import mpmath as mp
W = json.load(open(sys.argv[1])); mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 70; t0 = time.time()
L = mp.mpf(Fraction(W['L']).numerator) / Fraction(W['L']).denominator; k = int(W['EDGEK']); p = W['p']
cs = {int(j): Fraction(v) for j, v in W['coeffs'].items()}
D = max(cs) + 2 * k
coef = [mp.mpf(0)] * (D + 1)  # monomials in u = x/L
for j, c in cs.items():
    lp = mp.taylor(lambda u: mp.legendre(j, u) * (1 - u * u)**k, 0, j + 2 * k)
    for i, a in enumerate(lp): coef[i] += mp.mpf(c.numerator) / c.denominator * a
fu = lambda u: mp.polyval(coef[::-1], u)
f = lambda x: fu(x / L) if abs(x) <= L else mp.mpf(0)
# Gauss-Legendre nodes, Newton-refined
NG = D + 8
import numpy as np
xs, ws = np.polynomial.legendre.leggauss(NG); nodes = []
for x0 in xs:
    x = mp.mpf(x0)
    for _ in range(8):
        P, dP = mp.legendre(NG, x), mp.diff(lambda t: mp.legendre(NG, t), x); x -= P / dP
    w = 2 / ((1 - x * x) * mp.diff(lambda t: mp.legendre(NG, t), x)**2); nodes.append((x, w))
def g(x):
    a, b = x - L, L
    if a >= b: return mp.mpf(0)
    m, r = (a + b) / 2, (b - a) / 2
    return r * mp.fsum(w * f(m + r * t) * f(m + r * t - x) for t, w in nodes)
g0 = g(mp.mpf(0)); twoL = 2 * L
pole = 2 * mp.quad(lambda x: g(x) * 2 * mp.cosh(x / 2), [0, L, twoL])
a1 = mp.quad(lambda x: (g0 - g(x)) * mp.exp(-x) / mp.sinh(x), [0, L, twoL])
a2 = mp.quad(lambda x: g(x) * mp.exp(-x / 4) * mp.sinh(3 * x / 4) / mp.sinh(x), [0, L, twoL])
arch = (-mp.euler - mp.log(mp.pi)) * g0 - g0 * mp.log(1 - mp.exp(-4 * L)) + a1 - 2 * a2
def lam(n):
    for q in range(2, n + 1):
        if n % q == 0:
            m = n
            while m % q == 0: m //= q
            return q if m == 1 else None
prim = {}; n = 2
while mp.log(n) < twoL:
    q = lam(n)
    if q: prim[n] = 2 * mp.log(q) / mp.sqrt(n) * g(mp.log(n))
    n += 1
Qz = pole + arch - mp.fsum(prim.values())
Qko = Qz + mp.fsum(v for n, v in prim.items() if lam(n) == p)
print(json.dumps(dict(file=sys.argv[1], dps=mp.mp.dps, norm2=mp.nstr(g0, 15), pole=mp.nstr(pole, 20), arch=mp.nstr(arch, 20),
      Q_zeta=mp.nstr(Qz, 15), Q_ko=mp.nstr(Qko, 15), Q_ko_over_norm=mp.nstr(Qko / g0, 12), negative=bool(Qko < 0),
      wall_s=round(time.time() - t0, 1))), flush=True)
# Grok S57 (2026-10-17). RH OPEN.
