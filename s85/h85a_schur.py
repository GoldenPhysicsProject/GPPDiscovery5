"""Grok S85 H85a probe (NUMERICAL, float64, reuses s82os/os_hankel.build). Pole-sewing Schur scalar per parity.
Even: Q_e = R_e + p p^T with R_e (pole-free) of inertia (N-1,1)  =>  Q_e >= 0  <=>  s_e := 1 + p^T R_e^{-1} p <= 0  (PROVED, finite LA).
Odd:  Q_o = R_o - q q^T with R_o > 0                                =>  Q_o >= 0  <=>  s_o := 1 - q^T R_o^{-1} q >= 0.
When s = 0 the null vector is R^{-1}p (resp. R^{-1}q): the 'pole boundary state transferred by the pole-free resolvent'.
Test: overlap of R^{-1}p with the ground vector of Q. RH OPEN."""
import sys, json, math
sys.path.insert(0, "../s82os")
import numpy as np
from os_hankel import build
res = []
for L in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7):
    N = 28
    T, B = build("ZETA", L, N); c = build.comp
    p = c['ap'] + c['am']; q = c['ap'] - c['am']
    Qe = T + B; Qo = T - B
    Re_ = Qe - np.outer(p, p); Ro = Qo + np.outer(q, q)
    ere = np.linalg.eigvalsh(Re_); ero = np.linalg.eigvalsh(Ro)
    ye = np.linalg.solve(Re_, p); yo = np.linalg.solve(Ro, q)
    se = 1 + p @ ye; so = 1 - q @ yo
    eQe, VQe = np.linalg.eigh(Qe); eQo, VQo = np.linalg.eigh(Qo)
    ove = abs(ye @ VQe[:, 0]) / np.linalg.norm(ye); ovo = abs(yo @ VQo[:, 0]) / np.linalg.norm(yo)
    r = {"L": L, "N": N, "R_even_neg": int((ere < 0).sum()), "R_odd_neg": int((ero < 0).sum()),
         "s_even=1+p'R^-1p (<=0 iff even Weil)": se, "lam_min_even": eQe[0],
         "s_odd=1-q'R^-1q (>=0 iff odd Weil)": so, "lam_min_odd": eQo[0],
         "overlap(R^-1p, ground_even)": ove, "overlap(R^-1q, ground_odd)": ovo,
         "s_even*|R^-1p|^-2 vs lam_min (first order: lam ~ -s/|y|^2)": -se / (ye @ ye)}
    res.append(r); print(json.dumps(r), flush=True)
json.dump(res, open("h85a_schur_results.json", "w"), indent=1)
