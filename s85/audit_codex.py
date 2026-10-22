"""Grok S85: independent audit of Codex 6fc3833 (all-regulator boundary instability exponent) and 846ccf8
(critical GCD Mobius Euler packet). Exact rationals (fractions), sympy, mpmath. RH OPEN."""
import json, math
from fractions import Fraction as Fr
import mpmath as mp, sympy as sp, numpy as np
from sympy import primerange
out = {}
mp.mp.dps = 30
# ---- A1: kappa = 2 + digamma(1/4) - log pi  vs  2 - 2 zeta'/zeta(1/2) (pole of L[J0] at z=0) ----
k1 = 2 + mp.digamma(0.25) - mp.log(mp.pi)
k2 = 2 - 2 * mp.zeta(0.5, derivative=1) / mp.zeta(0.5)
out["A1_kappa"] = [str(k1), str(k2), float(abs(k1 - k2))]
# ---- A2: Laplace identities, symbolic ----
z, e, t, v = sp.symbols('z e t v', positive=True)
Ae = (z + e) * (z + e + sp.Rational(1, 2)) / (z * (z + sp.Rational(1, 2)))
pf = 1 + e * (1 + 2 * e) / z + e * (1 - 2 * e) / (z + sp.Rational(1, 2))
out["A2_transfer_partial_fractions_zero"] = str(sp.simplify(Ae - pf))
# Laplace of C_e(t) = 2 int_0^t (e^{-v/2}-e^{-t/2}) e^{(1-e)v} dv + 2(1-e^{-t/2}) at numeric e, z
def LC(ee, zz):
    f = lambda tt: 2*((mp.e**((0.5-ee)*tt)-1)/(0.5-ee) - mp.e**(-tt/2)*(mp.e**((1-ee)*tt)-1)/(1-ee)) + 2*(1-mp.e**(-tt/2))
    return mp.quad(lambda tt: mp.e**(-zz*tt) * f(tt), [0, 20, 60, mp.inf])
for ee, zz in ((0.3, 1.7), (0.8, 1.1)):
    claim = (1/(zz+ee-0.5) + 1) / (zz*(zz+0.5))
    out[f"A2_LaplaceC_e{ee}_z{zz}"] = [float(LC(ee, zz)), float(claim)]
# residue transfer: A_e(z_rho) * (-1/((rho-1/2)rho)) == -1/(z_rho(z_rho+1/2)), z_rho = rho-1/2-e
r = sp.symbols('r')
zr = r - sp.Rational(1, 2) - e
lhs = Ae.subs(z, zr) * (-1 / ((r - sp.Rational(1, 2)) * r)); rhs = -1 / (zr * (zr + sp.Rational(1, 2)))
out["A2_residue_transfer_zero"] = str(sp.simplify(lhs - rhs))
# ---- A3: explicit formula for J0 numerically (primes vs zeros) ----
def J0_primes(T):
    X = math.exp(T); ps = []
    from sympy import primerange as pr
    # J0 = int_1^X psi x^{-3/2} - int_1^X x^{-1/2} = sum_{n<=X} Lambda(n)*2(n^{-1/2}-X^{-1/2}) - 2(X^{1/2}-1)
    s = 0.0
    for p in pr(2, int(X) + 1):
        pk = p
        while pk <= X:
            s += math.log(p) * 2 * (pk**-0.5 - X**-0.5); pk *= p
    return s - 2 * (math.sqrt(X) - 1)
NZ = 400
mp.mp.dps = 15
gam = [float(mp.im(mp.zetazero(k))) for k in range(1, NZ + 1)]
mp.mp.dps = 30
def J0_zeros(T):
    kap = float(k1)
    R = 2*math.log(2*math.pi)*math.exp(-T/2) - sum(math.exp(-(2*n+0.5)*T)/(n*(4*n+1)) for n in range(1, 200))
    zs = 0.0
    for g in gam:
        rho = complex(0.5, g); zs += 2 * (np.exp((rho-0.5)*T) / (rho*(rho-0.5))).real
    return kap + R - zs
out["A3_J0_primes_vs_zeros(400 zeros; tail ~ sum 2/g^2 beyond)"] = {T: [J0_primes(T), J0_zeros(T)] for T in (3.0, 7.0, 11.0)}
out["A3_tail_bound_sum_2/g^2_beyond_g400"] = 2 * (math.log(gam[-1] / (2*math.pi)) + 1) / (2*math.pi*gam[-1])
# ---- A4: toy energy primitive + lower coefficient ----
d, g, T = sp.symbols('d g T', positive=True)
E = (sp.exp(2*d*T)-1)/(4*d) + sp.Rational(1, 2)*sp.re(sp.expand_complex((sp.exp((2*d+2*sp.I*g)*T)-1)/(2*d+2*sp.I*g)))
dd, gg, TT = 0.13, 2.7, 5.3
num = mp.quad(lambda x: mp.e**(2*dd*x)*mp.cos(gg*x)**2, [0, TT])
out["A4_toy_energy"] = [float(num), float(E.subs({d: dd, g: gg, T: TT}))]
# ---- B1: Mobius Euler packet G-norm exact = prod(1-1/p); optimality = 1/(G^{-1})_{11} on divisors of P# ----
from math import gcd
from itertools import combinations
B1 = []
P = list(primerange(2, 20))
for k in range(1, 8):
    Pk = P[:k]; divs = [1]
    for p in Pk: divs = divs + [d*p for d in divs]
    mu = lambda d: (-1)**sum(1 for p in Pk if d % p == 0)
    s = sum(Fr(mu(a)*mu(b), a*b//gcd(a, b)) for a in divs for b in divs)
    prod = 1
    for p in Pk: prod *= Fr(p-1, p)
    # prefix energy exact: sum a_m a_n / max(m,n), a_d = mu(d)/sqrt d (irrational) -> mpmath 40 digits
    pre = mp.fsum(mu(a)*mu(b)/mp.sqrt(a*b)/max(a, b) for a in divs for b in divs)
    B1.append({"k": k, "G_norm_exact": str(s), "prod(1-1/p)": str(prod), "equal": s == prod, "prefix_energy": float(pre), ">=1/2": bool(pre >= 0.5)})
out["B1"] = B1
# ---- B2: is the packet optimal beyond squarefree divisors? min_{a_1=1, supp<=N} ||a||_G^2 = 1/(G_N^{-1})_{11} ----
B2 = []
for N in (30, 100, 300, 1000, 2000):
    n = np.arange(1, N+1); Gm = np.gcd.outer(n, n) / np.sqrt(np.outer(n, n))
    e1 = np.zeros(N); e1[0] = 1
    Lc = np.linalg.cholesky(Gm); y = np.linalg.solve(Lc, e1); inv11 = y @ y
    pr = np.prod([1 - 1/p for p in primerange(2, N+1)])
    B2.append({"N": N, "min_G_norm_a1=1": 1/inv11, "prod_{p<=N}(1-1/p)": pr, "exp(-gamma)/logN": math.exp(-0.5772156649)/math.log(N)})
out["B2_min_over_n<=N"] = B2
# ---- B3: shifted products, first 200 primes ----
P200 = list(primerange(2, 1300))[:200]
out["B3_shift_200primes"] = {th: float(mp.fprod(1 - 2*mp.mpf(p)**(-1-th) + mp.mpf(p)**(-1-2*th) for p in P200)) for th in (-0.25, 0.0, 0.25)}
json.dump(out, open("audit_codex_results.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
