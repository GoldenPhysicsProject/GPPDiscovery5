"""Grok S88 H87a: does the Weil near-null vector align with the arch-only 'ghost' (negative) eigenvectors of A? High-precision NUMERICAL (mpmath), not interval. RH OPEN."""
import json, sys, mpmath as mp
from multiprocessing import Pool
sys.path.insert(0, "../s87"); import hp_forms as H
def legval(c, u):  # f(u)=sum c_i phi_i(u), phi_i orthonormal shifted Legendre on [0,1]
    return mp.fsum(c[i] * mp.sqrt(2 * i + 1) * mp.legendre(i, 2 * u - 1) for i in range(len(c)))
def nodes(c, M=400):
    v = [legval(c, mp.mpf(k) / M) for k in range(1, M)]
    return sum(1 for a, b in zip(v, v[1:]) if a * b < 0)
def task(args):
    FORM, th, L, N = args
    H.setup(90); d = H.build(FORM, L, N, th); out = {"form": FORM, "theta": th, "L": L, "N": N}
    for par in "eo":
        A = mp.matrix(d["A" + par]); P = mp.matrix(d["P" + par]); R = A - P; Q = R.copy()
        for p in d["p" + par]:
            pv = mp.matrix(p); Q += pv * pv.T
        EQ, VQ = mp.eigsy(Q); EA, VA = mp.eigsy(A)
        iq = min(range(N), key=lambda i: EQ[i]); v = VQ[:, iq]
        neg = sorted([i for i in range(N) if EA[i] < 0], key=lambda i: EA[i])
        ov = [float((v.T * VA[:, i])[0] ** 2) for i in neg]
        c = [v[i] for i in range(N)]
        o = {"qmin": mp.nstr(EQ[iq], 8), "A_neg": [mp.nstr(EA[i], 6) for i in neg], "ov_ghosts": ov, "ov_span": sum(ov),
             "vAv": mp.nstr((v.T * A * v)[0], 8), "vPv": mp.nstr((v.T * P * v)[0], 8), "nodes_v": nodes(c),
             "nodes_ghosts": [nodes([VA[j, i] for j in range(N)]) for i in neg]}
        if d["p" + par]:
            pv = mp.matrix(d["p" + par][0]); w = mp.lu_solve(R, pv); w = w / mp.norm(w)
            o["ov_Rinvp"] = float((v.T * w)[0] ** 2)
        out[par] = o
    return out
if __name__ == "__main__":
    Ls = [round(0.4 + 0.1 * i, 1) for i in range(11)]
    tasks = [("ZETA", 0, L, 24) for L in Ls] + [("DH", 0, L, 24) for L in Ls] + [("FT", 0.1, L, 24) for L in Ls]
    res = []
    with Pool(6) as pool:
        for r in pool.imap_unordered(task, tasks):
            res.append(r)
            print(r["form"], r["theta"], r["L"], "| even qmin", r["e"]["qmin"], "span", round(r["e"]["ov_span"], 4), "ov", [round(x, 4) for x in r["e"]["ov_ghosts"]], "nodes", r["e"]["nodes_v"], r["e"]["nodes_ghosts"],
                  "| odd qmin", r["o"]["qmin"], "span", round(r["o"]["ov_span"], 4), "ov", [round(x, 4) for x in r["o"]["ov_ghosts"]], "nodes", r["o"]["nodes_v"], r["o"]["nodes_ghosts"], flush=True)
    res.sort(key=lambda r: (r["form"], r["theta"], r["L"]))
    json.dump(res, open("ghost.json", "w"), indent=1)
