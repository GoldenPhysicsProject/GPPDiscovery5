"""Independent float reproduction of window negativity for F_theta (theta from the file), by summing
over ACTUAL zeta zeros (no explicit formula). NUMERICAL (float; zero sum truncated at N zeros).

Phi_g(s) = G(s-1/2) G(1/2-s), G(z) = int g e^{zx}; zeros of F_theta are rho +/- theta.
g(x) = sum_j c_j P_j(x/h) (1-(x/h)^2)^4 with rational c_j read from a coefficient file
(format: header 'h=<p/q> theta=<p/q> ...' then lines 'j c_j').

The coefficient files are Claude Code's (GPP-bridge, private, work/claude-code/2026-10-13_weil_window_certificate/)
and are NOT copied here; pass the path. The certified (Arb) values are Claude Code's; this script only
reproduces them in floating point.
Run: python experiments/ftheta_window/ftheta_zero_sum.py <coeff_file> [N=150]
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from fractions import Fraction
from mpmath import mp, mpf, quad, exp, legendre, zetazero, mpc, re
from common import write_result
mp.dps = 20
fn = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 150
L = open(fn).read().split("\n"); hdr = L[0].split()
h = Fraction(hdr[0][2:]); th = Fraction(hdr[1][6:])
c = {int(a): Fraction(b) for a, b in (l.split() for l in L[1:] if l.strip())}
hf, tf = mpf(h.numerator) / h.denominator, mpf(th.numerator) / th.denominator
g = lambda x: sum(mpf(v.numerator) / v.denominator * legendre(j, x / hf) for j, v in c.items() if v) * (1 - (x / hf) ** 2) ** 4
nodes = [-hf, -hf / 2, 0, hf / 2, hf]
G = lambda z: quad(lambda x: g(x) * exp(z * x), nodes)
norm = quad(lambda x: g(x) ** 2, nodes)
gam = [zetazero(n).imag for n in range(1, N + 1)]
def S(theta):
    tot = 0
    for gm in gam:
        for sh in ([theta, -theta] if theta else [0]):
            z = mpc(sh, gm); tot += 2 * re(G(z) * G(-z))
    return tot / norm
parity = "even" if all(j % 2 == 0 for j, v in c.items() if v) else ("odd" if all(j % 2 for j, v in c.items() if v) else "mixed")
Ft, Z = S(tf), S(0)
print(f"h={h} theta={th} N={N} parity={parity} F_theta={Ft} zeta={Z}")
write_result(os.path.join(os.path.dirname(__file__), f"result_h{h.numerator}_{h.denominator}.json"),
             {"experiment": "ftheta_window", "label": "NUMERICAL", "h": str(h), "theta": str(th), "N_zeros": N,
              "g_parity": parity, "Ftheta_sum_over_norm": str(Ft), "zeta_sum_over_norm": str(Z)},
             kind="float", method="direct sum over first N zeta zeros, mpmath quad")
