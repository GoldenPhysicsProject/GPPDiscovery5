"""Grok S97 run2: h-convergence (h=0.01) + drop-p locality controls (p=2,3,5,7). dps 60. RH OPEN."""
import sys, json
from multiprocessing import Pool
import run as R
if __name__ == "__main__":
    h = float(sys.argv[1])
    jobs = [("ZETA", "ZETA", 0, (), h, 1), ("FT0.1", "FT", 0.1, (), h, 1), ("DH", "DH", 0, (), h, 1)] + \
           [(f"ZETA-drop{p}", "ZETA", 0, (p,), h, 1) for p in (2, 3, 5, 7)]
    with Pool(7) as p: out = p.map(R.task, jobs)
    json.dump(out, open(f"chain2_h{h}.json", "w"), indent=0)
    for r in out: print(r["name"], "n", r["n"], "fail_k", r["fail_k"], "T_fail", r["T_fail"], {k: r[k] for k in r if k.startswith("chol")})
