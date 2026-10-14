"""Li coefficients lambda_n for xi, F_theta and the completed Davenport-Heilbronn function.
NUMERICAL (float, mpmath; Cauchy-integral Taylor coefficients, no rigorous error bound).

Definition (Li 1997; Bombieri-Lagarias 1999): for an entire F with F(s)=F(1-s), F(1) != 0,
  lambda_n = sum_rho [1 - (1 - 1/rho)^n]   (zeros with multiplicity, symmetric summation),
equivalently  log F(1/(1-z)) = log F(1) + sum_{n>=1} (lambda_n / n) z^n  near z = 0.
We get lambda_n = n * [z^n] log F(1/(1-z)) by the trapezoid rule on |z| = r (two radii as a check).
The Taylor radius is 1 / max|1 - 1/rho|; it is 1 under RH and < 1 when F has zeros with Re rho < 1/2
(F_theta, DH), so r must stay below it. Normalisation: F = xi (xi(1)=1/2); F_theta = xi(s-theta)xi(s+theta)
(so F_0 = xi^2 has lambda_n doubled); DH: Lambda(s) = (5/pi)^{s/2} Gamma((s+1)/2) f(s).
Run: python experiments/li_coefficients/li_coeffs.py [nmax=100]
"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from mpmath import mp, mpf, mpc, exp, log, pi, nint
from common import xi, dh_completed, write_result
mp.dps = 60
nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 100
M = 768
def taylor_log(F, r):
    vals = []
    prev = None
    for k in range(M):
        z = r * exp(2j * pi * k / M)
        v = log(F(1 / (1 - z)))
        if prev is not None:  # continuous branch along the circle
            v += 2j * pi * nint((prev - v).imag / (2 * pi))
        vals.append(v); prev = v
    assert abs((vals[-1] - vals[0]).imag) < 1, "winding: a zero lies inside |z|<r"
    return [sum(vals[k] * exp(-2j * pi * k * n / M) for k in range(M)) / M / r**n for n in range(nmax + 1)]
funcs = {
    "xi": lambda s: xi(s),
    "F_theta=1/4": lambda s: xi(s - mpf(1)/4) * xi(s + mpf(1)/4),
    "F_theta=1/20": lambda s: xi(s - mpf(1)/20) * xi(s + mpf(1)/20),
    "DH": dh_completed,
}
which = sys.argv[2:] or list(funcs)
out = {}
for name in which:
    F = funcs[name]
    a, b = taylor_log(F, mpf("0.5")), taylor_log(F, mpf("0.4"))
    lam = [n * a[n].real for n in range(1, nmax + 1)]
    dif = max(abs(n * (a[n] - b[n])) for n in range(1, nmax + 1))
    imax = max(abs(n * a[n].imag) for n in range(1, nmax + 1))
    neg = [n for n, l in enumerate(lam, 1) if l < 0]
    out[name] = {"lambda": [float(x) for x in lam], "radius_check_maxdiff": float(dif),
                 "max_imag_part": float(imax), "negative_n": neg}
    print(name, "lambda_1..5", [round(float(x), 6) for x in lam[:5]], "lambda_nmax", float(lam[-1]),
          "neg:", neg[:10], "r-check", float(dif), flush=True)
tag = "_".join(w.replace("=", "").replace("/", "_") for w in which) if sys.argv[2:] else "all"
write_result(os.path.join(os.path.dirname(__file__), f"result_{tag}.json"),
             {"experiment": "li_coefficients", "label": "NUMERICAL", "nmax": nmax, "M": M, "radii": [0.5, 0.4], "data": out},
             kind="float", method="trapezoid Cauchy integral of log F(1/(1-z)), mpmath dps 60")
