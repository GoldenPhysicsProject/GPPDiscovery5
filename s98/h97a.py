"""Grok S98 -- H97a test: sign of each prime-power atom of the Weil window form on the soft mode.
Uses the S97 Toeplitz Gram data (s97/chain2_h{h}.json, dps 60 entries, 20 digits stored) for ZETA, DH and ZETA-drop-p.
Soft mode v = lowest eigenvector of the step-basis Toeplitz Gram on [0,T] (float64 eigh; cond <= 7e4, inputs 20 digits).
Atom n contributes -2 a_n n^{-1/2} g_v(log n), g_v = autocorrelation of the step function v (tent kernel).
Also Q_p(v) = v^T (T_zeta - T_drop_p) v near T = log p^+. NUMERICAL, not interval. RH OPEN."""
import json, math, sys, numpy as np
sys.path.insert(0, "../s87"); import hp_forms as H
from scipy.linalg import toeplitz, eigh
out = []
for h in (0.02, 0.01):
    R = {r["name"]: r for r in json.load(open(f"../s97/chain2_h{h}.json"))}
    cz = np.array([float(x) for x in R["ZETA"]["c"]])
    out.append(f"== h={h}: Q_p on zeta soft mode just past log p ==")
    for p in (2, 3, 5, 7):
        cd = np.array([float(x) for x in R[f"ZETA-drop{p}"]["c"]]); cp = cz - cd
        for d in (0.02, 0.06, 0.2, 0.4):
            n = int(round((math.log(p) + d) / h)); w, V = eigh(toeplitz(cz[:n])); v = V[:, 0]
            wd = eigh(toeplitz(cd[:n]), eigvals_only=True)[0]
            nodeless = bool(np.all(v > 0) or np.all(v < 0))
            out.append(f" p={p} T={n*h:.2f} lam_min={w[0]:.3e} Q_p(v)={v@toeplitz(cp[:n])@v:+.3e} "
                       f"Q_p/lam={v@toeplitz(cp[:n])@v/w[0]:+.2f} nodeless={nodeless} drop_lam_min={wd:+.3e}")
    dw = H.dh_weights(30)
    def g(v, z):
        n = len(v); r = np.array([v[:n - m] @ v[m:] for m in range(n)]); m = np.arange(n)
        return float(np.sum(r * np.maximum(h - np.abs(z - m * h), 0)))
    out.append(f"== h={h}: per-atom contribution -2 a_n g_v(log n)/sqrt(n) on the soft mode ==")
    for nm in ("ZETA", "DH"):
        c = np.array([float(x) for x in R[nm]["c"]])
        for T in (1.5, 2.2, 3.0):
            n = int(round(T / h)); w, V = eigh(toeplitz(c[:n])); v = V[:, 0]
            nodeless = bool(np.all(v > 0) or np.all(v < 0)); rows = []
            for k in range(2, int(math.exp(T)) + 1):
                a = float(H.mang(k)) if nm == "ZETA" else float(dw.get(k, 0))
                if a == 0: continue
                gk = g(v, math.log(k)); rows.append(f"{k}:a={a:+.3f},g{'+' if gk > 0 else '-'},{-2*a/math.sqrt(k)*gk:+.1e}")
            out.append(f" {nm} T={T} lam_min={w[0]:.2e} nodeless={nodeless} " + " ".join(rows))
open("h97a.out", "w").write("\n".join(out) + "\n"); print("\n".join(out))
