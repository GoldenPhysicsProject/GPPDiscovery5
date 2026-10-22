"""Grok S87 H86b scan: inertia of arch-only A and pole-free R=A-P, both parities, ZETA/DH/F_theta, high precision."""
import json, sys, mpmath as mp
from multiprocessing import Pool
import hp_forms as H
def task(args):
    FORM, th, L, N = args
    H.setup(90 if N <= 24 else 110)
    d = H.build(FORM, L, N, th); out = {"form": FORM, "theta": th, "L": L, "N": N}
    for par in "eo":
        A = d["A" + par]; R = H.add(A, d["P" + par], 1, -1)
        eA, eR = H.eigs(A), H.eigs(R)
        out[par] = {"nA": sum(1 for x in eA if x < 0), "nR": sum(1 for x in eR if x < 0),
                    "A_low": [mp.nstr(x, 8) for x in eA[:5]], "R_low": [mp.nstr(x, 8) for x in eR[:4]]}
    return out
if __name__ == "__main__":
    Ls = [round(0.30 + 0.05 * i, 2) for i in range(23)]
    tasks = [("ZETA", 0, L, N) for L in Ls for N in (24, 32)] + [(f, th, L, 24) for L in Ls for f, th in (("DH", 0), ("FT", 0.1), ("FT", 0.25))]
    tasks.sort(key=lambda t: -t[3])
    res = []
    with Pool(8) as pool:
        for r in pool.imap_unordered(task, tasks):
            res.append(r); print(r["form"], r["theta"], r["L"], r["N"], "even nA", r["e"]["nA"], "nR", r["e"]["nR"], "| odd nA", r["o"]["nA"], "nR", r["o"]["nR"], flush=True)
    res.sort(key=lambda r: (r["form"], r["theta"], r["N"], r["L"]))
    json.dump(res, open("scan.json", "w"), indent=1)
