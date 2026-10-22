"""Grok S91: cache the S87 high-precision window forms (mpmath dps 90) for s91 (h90a.py, schur_audit.py).
Per (FORM, theta, L, N): A_e, A_o (arch), P_e, P_o (prime sums), pole vectors p (even), q (odd), stored as mp strings.
Conventions (CHECKED s89/calib.py): Q_e = A_e - P_e + p p^T,  Q_o = A_o - P_o - q q^T.  HIGH-PRECISION NUMERICAL, not interval. RH OPEN."""
import sys, os, pickle, mpmath as mp
from multiprocessing import Pool
sys.path.insert(0, "../s87"); import hp_forms as H
def task(a):
    F, th, L, N = a
    H.setup(90); d = H.build(F, L, N, th)
    st = lambda M: [[mp.nstr(x, 95) for x in r] for r in M]
    return (F, th, L, N), {k: st(d[k]) for k in ("Ae", "Ao", "Pe", "Po")} | {"pe": [[mp.nstr(x, 95) for x in v] for v in d["pe"]], "po": [[mp.nstr(x, 95) for x in v] for v in d["po"]]}
def load(F, th, L, N, dps=90):
    mp.mp.dps = dps
    cache = pickle.load(open("forms_cache.pkl", "rb"))
    c = cache[(F, th, L, N)] if (F, th, L, N) in cache else task((F, th, L, N))[1]
    mp.mp.dps = dps
    M = lambda k: mp.matrix([[mp.mpf(x) for x in r] for r in c[k]])
    out = {k: M(k) for k in ("Ae", "Ao", "Pe", "Po")}
    out["pe"] = [mp.matrix([mp.mpf(x) for x in v]) for v in c["pe"]]; out["po"] = [mp.matrix([mp.mpf(x) for x in v]) for v in c["po"]]
    return out
def Qform(d, par):
    Q = d["A" + par] - d["P" + par]; sg = 1 if par == "e" else -1
    for v in d["p" + par]: Q += sg * v * v.T
    return Q
if __name__ == "__main__":
    LA = [0.6, 0.8, 1.0, 1.2]; LB = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    tasks = [("ZETA", 0, L, 24) for L in LB] + [("DH", 0, L, 24) for L in LA] + [("FT", 0.1, L, 24) for L in LA] \
          + [("ZETA", 0, L, 32) for L in (0.6, 1.0, 1.2, 1.4)]
    with Pool(6) as pool: res = dict(pool.map(task, tasks))
    pickle.dump(res, open("forms_cache.pkl", "wb")); print("cached", len(res))
