"""Grok S99 -- H98a: copositivity of window forms on the cone v >= 0 (step basis, S97 dps-60 Toeplitz data, 20 digits).
For each form and window T = n*h: lam_min (float64 eigh) vs copositive minimum mu_+ = min_{v>=0,|v|=1} v^T T v
(projected gradient + 40 random restarts + eigenvector-positive-part starts). NUMERICAL, not interval. RH OPEN."""
import json, numpy as np
from scipy.linalg import toeplitz, eigh
rng = np.random.default_rng(99)
def copos(A, iters=4000):
    w, V = eigh(A); L = w[-1]; best = np.inf
    starts = [np.abs(V[:, 0]), np.maximum(V[:, 0], 0), np.maximum(-V[:, 0], 0), np.ones(len(A))] + [rng.random(len(A)) for _ in range(40)]
    for v in starts:
        v = v / (np.linalg.norm(v) + 1e-300)
        for _ in range(iters):
            g = A @ v; v2 = np.maximum(v - g / L, 0); nv = np.linalg.norm(v2)
            if nv == 0: break
            v2 /= nv
            if np.linalg.norm(v2 - v) < 1e-13: v = v2; break
            v = v2
        q = v @ A @ v
        if q < best: best, bv = q, v
    return best, bv
out = []
for h in (0.02,):
    R = {r["name"]: r for r in json.load(open(f"../s97/chain2_h{h}.json"))}
    for nm, r in R.items():
        c = np.array([float(x) for x in r["c"]])
        for n in (25, 40, 50, 60, 75, 90, 110, 130, 150):
            A = toeplitz(c[:n]); w, V = eigh(A); v0 = V[:, 0]
            sc = int(np.sum(np.diff(np.sign(v0[np.abs(v0) > 1e-12])) != 0))
            mu, bv = copos(A)
            supp = int(np.sum(bv > 1e-8))
            out.append(f"{nm:11s} T={n*h:.2f} lam_min={w[0]:+.3e} signchg(v0)={sc} copos_min={mu:+.3e} supp={supp}/{n} "
                       f"{'PD' if w[0] > 0 else 'NOT-PD'} {'COPOS' if mu > 0 else 'NOT-COPOS'}")
open("h98a.out", "w").write("\n".join(out) + "\n"); print("\n".join(out))
