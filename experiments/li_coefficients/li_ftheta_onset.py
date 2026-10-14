"""Where do the Li coefficients of F_theta first go negative? NUMERICAL (float) + asymptotic estimate.

lambda_n(F_theta) = 2 lambda_n(xi) + D_n, with the EXACT zero-side identity
  D_n = sum_{rho of zeta, Im>0, conj pairs} 2 Re[ 2 w0^n - w_+^n - w_-^n ],
  w0 = 1 - 1/rho,  w_pm = 1 - 1/(rho -/+ theta)   (zeros of F_theta are rho + theta and rho - theta).
D_n converges fast (terms ~ gamma^-4 for gamma^2 >> n theta), so we sum the first N zeta zeros.
2 lambda_n(xi) is taken from the asymptotic lambda_n(xi) = (n/2)(log n - 1 + gamma_E - log 2pi) + O(sqrt(n) log n)
(Voros 2006; Coffey). The O(sqrt n log n) error makes the onset estimate HEURISTIC, not certified.
Run: python experiments/li_coefficients/li_ftheta_onset.py [theta=0.25] [N=600]
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from mpmath import mp, mpf, mpc, zetazero, log, euler, pi, re
from common import write_result
mp.dps = 30
theta = mpf(sys.argv[1]) if len(sys.argv) > 1 else mpf("0.25")
N = int(sys.argv[2]) if len(sys.argv) > 2 else 600
rhos = [zetazero(k) for k in range(1, N + 1)]
W = [(1 - 1 / r, 1 - 1 / (r - theta), 1 - 1 / (r + theta)) for r in rhos]
def D(n):
    return sum(2 * re(2 * w0**n - wp**n - wm**n) for w0, wp, wm in W)
def lam_xi_asym(n):
    return n / mpf(2) * (log(n) - 1 + euler - log(2 * pi))
rows = []; first = None
for n in list(range(500, 30001, 250)):
    d = D(n); tot = 2 * lam_xi_asym(n) + d
    rows.append((n, float(d), float(2 * lam_xi_asym(n)), float(tot)))
    if first is None and tot < 0: first = n
print("theta", theta, "N zeros", N, "first n (step 250) with 2*lam_asym(xi)+D_n < 0:", first)
for r in rows[::8]: print(r)
write_result(os.path.join(os.path.dirname(__file__), f"result_onset_theta{float(theta)}.json"),
             {"experiment": "li_ftheta_onset", "label": "NUMERICAL (heuristic: xi part asymptotic)",
              "theta": str(theta), "N_zeros": N, "first_negative_n_step250": first,
              "rows_n_D_2lamxi_total": rows}, kind="float", method="zero sums (mpmath zetazero) + Voros asymptotic")
