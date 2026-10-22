"""Grok S86 H86b test: inertia of archimedean-only A_L (even+odd) vs number of prime powers n<=e^{2L} (Hankel entries), zeta. NUMERICAL float64, N=32."""
import json, math, numpy as np
from h86_prime_path import build
from sympy import factorint
pp = lambda X: sum(1 for n in range(2, int(X) + 1) if len(factorint(n)) == 1)
rows = []
for L in np.arange(0.30, 1.251, 0.05):
    d = build("ZETA", float(L), 32, M=80, S=128)
    ne = int((np.linalg.eigvalsh(d["Ae"]) < 0).sum()); no = int((np.linalg.eigvalsh(d["Ao"]) < 0).sum())
    nRe = int((np.linalg.eigvalsh(d["Ae"] - d["Pe"]) < 0).sum()); nRo = int((np.linalg.eigvalsh(d["Ao"] - d["Po"]) < 0).sum())
    r = {"L": round(float(L), 3), "nA_e": ne, "nA_o": no, "nR_e": nRe, "nR_o": nRo, "#pp<=e^2L": pp(math.exp(2 * L)), "#pp<=e^L": pp(math.exp(L))}
    rows.append(r); print(json.dumps(r), flush=True)
json.dump(rows, open("h86_inertia_count.json", "w"), indent=1)
