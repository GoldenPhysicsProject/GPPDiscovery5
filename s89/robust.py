"""S89 robustness: H_t onset points with M=40 swapped zeros (vs 50) and N=24 (vs 32); zeta odd rho_out-1 digits."""
import os, json, sys
os.environ["MZ"] = "40"
import band_split as B, mpmath as mp
from multiprocessing import Pool
def t(a):
    r = B.task(a); return {"args": a, "M": r.get("M_zeros"), "e": {k: r["e"][k] for k in ("qmin", "out_min", "in_schur_min", "nodes_v", "span_v")},
                           "o": {k: r["o"][k] for k in ("qmin", "out_min", "in_schur_min", "nodes_v", "span_v")}}
if __name__ == "__main__":
    tasks = [("HT", "-2", 0.7, 32), ("HT", "-2", 0.8, 32), ("HT", "-4", 0.5, 32), ("HT", "-4", 0.6, 32), ("HT", "-2", 0.7, 24), ("HT", "-4", 0.5, 24)]
    with Pool(6) as p:
        out = p.map(t, tasks)
    for o in out: print(json.dumps(o))
    json.dump(out, open("robust.json", "w"), indent=1)
