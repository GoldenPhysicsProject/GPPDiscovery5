"""Grok S86: N-robustness of the inertia of A (arch only) vs R = A - P (pole-free) for zeta. NUMERICAL float64."""
import json, numpy as np
from h86_prime_path import build
res = []
for L in (0.35, 0.4, 0.5, 0.7):
    for N in (16, 24, 32):
        d = build("ZETA", L, N); r = {"L": L, "N": N}
        for par in ("e", "o"):
            A, P = d["A" + par], d["P" + par]; eA = np.linalg.eigvalsh(A); eR = np.linalg.eigvalsh(A - P)
            r[par] = {"A_low3": eA[:3].tolist(), "R_low2": eR[:2].tolist(), "nA": int((eA < 0).sum()), "nR": int((eR < 0).sum())}
        res.append(r); print(json.dumps(r), flush=True)
json.dump(res, open("h86_inertia_N.json", "w"), indent=1)
