"""Parity split of the Weil pole term. NUMERICAL (float quadrature, mpmath).

For real g supported in [-1,1], f = g * g~ (autocorrelation), fhat(s) = int f(x) e^{(s-1/2)x} dx,
G(z) = int g(x) e^{zx} dx. Claim checked: fhat(0) + fhat(1) = 2 G(1/2) G(-1/2)
= 2(C^2 - S^2) with C = int g cosh(x/2), S = int g sinh(x/2); hence +2G(1/2)^2 for even g,
-2G(1/2)^2 for odd g, and no square for mixed parity.
Run: python experiments/pole_parity/pole_parity.py   (writes result.json next to the script)
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from mpmath import mp, mpf, quad, exp, cosh, sinh
from common import write_result
mp.dps = 25
h = 1
def run(g):
    G = lambda z: quad(lambda x: g(x) * exp(z * x), [-h, h])
    f = lambda x: quad(lambda y: g(y) * g(y - x), [max(-h, x - h), min(h, x + h)]) if abs(x) < 2 * h else 0
    fh = lambda s: quad(lambda x: f(x) * exp((s - mpf(1) / 2) * x), [-2 * h, 0, 2 * h])
    C = quad(lambda x: g(x) * cosh(x / 2), [-h, h]); S = quad(lambda x: g(x) * sinh(x / 2), [-h, h])
    return dict(pole=fh(0) + fh(1), two_G_half_G_minus_half=2 * G(0.5) * G(-0.5),
                two_C2_minus_S2=2 * (C**2 - S**2), two_G_half_sq=2 * G(0.5) ** 2)
cases = {"even (1-x^2)^3": lambda x: (1 - x * x) ** 3,
         "odd x(1-x^2)^3": lambda x: x * (1 - x * x) ** 3,
         "mixed (1-x^2)^3(1+2x)": lambda x: (1 - x * x) ** 3 * (1 + 2 * x)}
res = {k: {a: str(b) for a, b in run(g).items()} for k, g in cases.items()}
for k, v in res.items(): print(k, v)
write_result(os.path.join(os.path.dirname(__file__), "result.json"),
             {"experiment": "pole_parity", "label": "NUMERICAL", "cases": res},
             kind="float", method="mpmath tanh-sinh quadrature, nested")
