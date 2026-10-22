import sys, time, numpy as np
sys.path.insert(0, "../s86"); sys.path.insert(0, "../s82os")
import hp_forms as H
from h86_prime_path import build as fb
H.setup(60)
for FORM, th in (("ZETA", 0), ("DH", 0), ("FT", 0.25)):
    t = time.time(); d = H.build(FORM, 0.6, 12, th); t1 = time.time() - t
    f = fb(FORM, 0.6, 12, th)
    for k in ("Ae", "Ao", "Pe", "Po"):
        print(FORM, k, "maxdiff", np.max(np.abs(np.array(d[k], dtype=float) - f[k])), "time", round(t1, 2))
    if d["pe"]: print("pe diff", np.max(np.abs(np.array(d["pe"], float).T - f["pe"])))
