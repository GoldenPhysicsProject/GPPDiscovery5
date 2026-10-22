"""Grok S86: independent audit of Codex 0f238bc/38c633a (GCD/Poisson coefficient recovery threshold).
K_eps(m,n) = (gcd(m,n)^2/(mn))^eps on positive integers.
Claims: (i) ||coeff_n||^2_{H*} = zeta(2eps) prod_{p|n}(1+p^{-2eps}) (eps>1/2);
(ii) ||B_T||^2 = (1-e^{-T}) zeta(2eps) for 0<T<=log2; (iii) eps<=1/2: Mobius packet norm prod(1-p^{-2eps}) -> 0, image fixed.
Method: exact Fractions on tensor boxes (eps=1/2 gives rationals after squaring? no: use eps with p^{-2eps} rational => eps=1/2 or 1),
float64 for general eps on divisor-closed sets {n<=N}. NUMERICAL/EXACT as labelled."""
import math, json, itertools
from fractions import Fraction
import numpy as np
from sympy import primerange, factorint
from mpmath import zeta
out = {}
def gram(ns, eps):
    ns = np.array(ns, dtype=float)
    G = np.array([[ (math.gcd(int(a), int(b))**2/(a*b))**eps for b in ns] for a in ns])
    return G
# (i) exact at eps=1 (r^2 = p^-2 rational): tensor box primes {2,3,5}, exponents<=3 ; dual norm = (K^-1)_{nn}
def exact_box(primes, emax, eps2):  # eps2 = 2*eps integer -> K entries rational p^{-eps|i-j|}: need eps integer -> eps=1
    ns = [math.prod(p**e for p, e in zip(primes, es)) for es in itertools.product(range(emax+1), repeat=len(primes))]
    import sympy
    K = sympy.Matrix(len(ns), len(ns), lambda i, j: sympy.Rational(math.gcd(ns[i], ns[j])**2, ns[i]*ns[j]))  # eps=1
    Ki = K.inv()
    return ns, Ki
ns, Ki = exact_box([2, 3], 3, 2)
chk = []
for n in (1, 2, 3, 6, 4, 12, 8):
    i = ns.index(n)
    fac = factorint(n); pred = 1
    for p in (2, 3):
        r2 = Fraction(1, p*p)
        e = fac.get(p, 0)
        pred *= (1/(1-r2)) if (e == 0 or e == 3) else (1+r2)/(1-r2)  # e==emax is terminal: (K^-1) corner = 1/(1-r^2)
    chk.append({"n": n, "Kinv_nn": str(Ki[i, i]), "pred_box": str(pred), "equal": Fraction(str(Ki[i, i])) == pred})
out["exact_eps1_box_{2,3}^3"] = chk
# (i)/(ii) float: coefficient-1 dual norm over all n<=N vs zeta(2eps); eps in {0.75, 0.6, 0.5}
rows = []
for eps in (0.75, 0.6, 0.5):
    for N in (64, 256, 1024):
        G = gram(range(1, N+1), eps)
        e1 = np.zeros(N); e1[0] = 1
        v = np.linalg.solve(G, e1)
        Pn = math.prod(1/(1-p**(-2*eps)) for p in primerange(2, N+1))
        rows.append({"eps": eps, "N": N, "(K^-1)_11": float(v[0]), "Euler_prod_p<=N": Pn,
                     "zeta(2eps)": float(zeta(2*eps)) if eps > 0.5 else "inf", "cond": float(np.linalg.cond(G))})
out["coeff1_dual_norm"] = rows
# (ii) B_T norm on T<=log2: B_T F = e^{-t/2} a_1 => ||B_T||^2 = (1-e^{-T}) * ||coeff_1||^2 (trivial given (i))
# (iii) Mobius packet at eps=1/2 exact source norm on primorial divisors
pk = []
for k in (3, 5, 7):
    P = list(primerange(2, 100))[:k]
    ds = [math.prod(c) for r in range(k+1) for c in itertools.combinations(P, r)]
    a = {d: (-1)**len(factorint(d)) * d**-0.5 for d in ds}
    G = gram(ds, 0.5); av = np.array([a[d] for d in ds])
    pk.append({"k": k, "norm2": float(av @ G @ av), "prod(1-1/p)": math.prod(1-1/p for p in P)})
out["mobius_packet_eps_half"] = pk
json.dump(out, open("audit_codex_threshold.json", "w"), indent=1, default=str); print(json.dumps(out, indent=1, default=str))
