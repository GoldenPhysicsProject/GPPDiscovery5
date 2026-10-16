#!/usr/bin/env python3
"""H40a.B: Spot-check S(t) assembly and top-hat identity caveats.

S(t)=W(t)-A(t) with W=Pole_t2+Arch_t2 using exact sum 1/a_j^2 = pi^2/4+2G.
Identity S=(t/2)Q_zeta(f_t): CONFIRMED-with-caveat (mollification).
"""
from __future__ import annotations
import json, os, math
import numpy as np
from mpmath import mp, mpf, nstr, cosh, exp, log, euler, pi, catalan

OUT = os.path.join(os.path.dirname(__file__), "h40a_tophat.json")
mp.dps = 50
omega0 = -euler - pi/2 - 3*log(2) - log(pi)
sum_inv_a2 = pi**2/4 + 2*catalan

def assemble(t, J=200):
    t = mpf(t)
    Nmax = max(2, int(math.floor(math.exp(float(t)))) + 2)
    sieve = np.ones(Nmax + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(math.isqrt(Nmax)) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    lam = np.zeros(Nmax + 1, dtype=float)
    for p in np.flatnonzero(sieve):
        n = int(p)
        while n <= Nmax:
            lam[n] = math.log(int(p))
            n *= int(p)
    A = mpf(0)
    for n in range(2, Nmax + 1):
        if lam[n] == 0:
            continue
        un = log(n)
        if un < t:
            A += mpf(lam[n]) / mp.sqrt(n) * (t - un)
    Pole = 8 * (cosh(t / 2) - 1)
    Arch_lin = (omega0 / 2) * t
    sum_exp = mpf(0)
    for j in range(J):
        a = 2 * j + mpf("0.5")
        sum_exp += exp(-a * t) / (a * a)
    aJ = 2 * J + mpf("0.5")
    z4 = exp(-2 * t)
    r = min(z4, mpf("0.99"))
    tail = (exp(-aJ * t) / (aJ * aJ)) / (1 - r)
    # gamma = sum(1-e)/a^2 ∈ [sum1 - sum_exp - tail, sum1 - sum_exp]
    gamma_lo = sum_inv_a2 - sum_exp - tail
    gamma_hi = sum_inv_a2 - sum_exp
    Arch_lo = Arch_lin + gamma_lo
    Arch_hi = Arch_lin + gamma_hi
    W_lo, W_hi = Pole + Arch_lo, Pole + Arch_hi
    S_lo, S_hi = W_lo - A, W_hi - A
    return {
        "t": nstr(t, 20),
        "Pole_t2": nstr(Pole, 25),
        "A_Prime_t2": nstr(A, 25),
        "S_enclosure": [nstr(S_lo, 25), nstr(S_hi, 25)],
        "S_mid": nstr((S_lo + S_hi) / 2, 25),
        "Q_proxy_S_over_t_half": nstr(((S_lo + S_hi) / 2) / (t / 2), 25),
        "tail": nstr(tail, 10),
    }

spots = [assemble(tv) for tv in [1.0, float(mp.log(2)), 2.0, 5.0]]
# At t=log2, expect S≈0.06355558603757 (Codex / our cone)
t0_ok = abs(float(spots[1]["S_mid"]) - 0.06355558603757) < 1e-12

results = {
    "spots": spots,
    "t0_matches_cone": t0_ok,
    "piece_defs": "CONFIRMED (W,A by definition)",
    "identity_S_eq_t_half_Q": {
        "status": "CONFIRMED-with-caveat",
        "note": (
            "By explicit formula, (t/2)Pole/Arch/Prime equal the displayed pieces "
            "for the top-hat autocorrelation g=(1-|u|/t)_+. Top-hat is L2 discontinuous; "
            "mollify and pass gamma pairing by dominated convergence "
            "(sinc^2(r) log(2+|r|) integrable); finite prime sum passes. "
            "GAP only under rigid Schwartz-only statements of Weil. "
            "No zero-sum numerical Q check (would not audit independently)."
        ),
    },
    "label": "CONFIRMED-with-caveat",
}
with open(OUT, "w") as f:
    json.dump(results, f, indent=2)
print("t0_ok", t0_ok)
for s in spots:
    print(s["t"], s["S_mid"])
print("Wrote", OUT)
