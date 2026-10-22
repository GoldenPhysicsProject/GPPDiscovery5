"""Grok S86: L at which the archimedean-only parity form A_L gains a negative eigenvalue (zeta, N=16, bisection). NUMERICAL float64.
Compare with prime-power entry thresholds: Hankel part n enters at L = log(n)/2, Toeplitz part at L = log(n)."""
import json, math, numpy as np
from h86_prime_path import build
def nA(par, L, k):
    d = build("ZETA", L, 16, M=64, S=96); return np.linalg.eigvalsh(d["A" + par])[k]
out = {}
for par, k, lo, hi in (("o", 0, 0.35, 0.4), ("e", 1, 0.5, 0.7)):
    for _ in range(18):
        mid = (lo + hi) / 2
        if nA(par, mid, k) < 0: hi = mid
        else: lo = mid
    out[f"{par}: eigenvalue #{k} of A crosses 0 at L"] = (lo + hi) / 2
out["prime entries (L)"] = {"Hankel n=2": math.log(2)/2, "Hankel n=3": math.log(3)/2, "Hankel n=4": math.log(4)/2, "Toeplitz n=2": math.log(2)}
print(json.dumps(out, indent=1)); json.dump(out, open("h86_arch_crossings.json", "w"), indent=1)
