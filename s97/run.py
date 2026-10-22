"""Grok S97 run: Levinson chain (Schur/Verblunsky kappa_k, pivots P_k) for ZETA, DH, FT0.1, FT0.05, ZETA-drop2, ZETA-flip-scale.
Cross-check PD by mp Cholesky of the full Toeplitz at the last PD n and LDL^T inertia at the end. dps 60. RH OPEN."""
import sys, json, mpmath as mp
from multiprocessing import Pool
import chain as C
T = 3.0
def task(a):
    name, F, th, drop, h, ps = a
    C.setup(60); S = C.spec(F, th, drop=drop); n = int(round(T / h))
    c = C.toeplitz_c(S, h, n, ps); kap, P, fail = C.levinson(c)
    # independent check: Cholesky on full T_n and on T_fail
    def chol(m):
        M = mp.matrix(m, m)
        for i in range(m):
            for j in range(m): M[i, j] = c[abs(i - j)]
        try: mp.cholesky(M); return True
        except Exception: return False
    chk = {"chol_full": chol(n)}
    if fail: chk.update({"chol_before_fail": chol(fail), "chol_at_fail": chol(fail + 1)})
    return {"name": name, "h": h, "n": n, "c": [mp.nstr(x, 20) for x in c], "kappa": [mp.nstr(x, 20) for x in kap],
            "P": [mp.nstr(x, 20) for x in P], "fail_k": fail, "T_fail": None if fail is None else (fail + 1) * h, **chk}
if __name__ == "__main__":
    h = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
    jobs = [("ZETA", "ZETA", 0, (), h, 1), ("DH", "DH", 0, (), h, 1), ("FT0.1", "FT", 0.1, (), h, 1), ("FT0.05", "FT", 0.05, (), h, 1),
            ("ZETA-drop2", "ZETA", 0, (2,), h, 1), ("ZETA-pole0.9", "ZETA", 0, (), h, 0.9)]
    with Pool(6) as p: out = p.map(task, jobs)
    json.dump(out, open(f"chain_h{h}.json", "w"), indent=0)
    for r in out: print(r["name"], "n", r["n"], "fail_k", r["fail_k"], "T_fail", r["T_fail"], {k: r[k] for k in r if k.startswith("chol")})
