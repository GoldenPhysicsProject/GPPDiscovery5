"""Grok S92 -- H_-t (backward heat flow) on the Connes-Consani near-radical span S = span E(prolates) (the S91 H90a span).
Q_HT(t) = Q_zeta + Delta_t (s89/band_split.delta_ht: zero-side swap of the first MZ zeros of zeta for those of H_t, t = -1,-2,-4).
Question: does zeta's near-radical span WITNESS H_-t's failure (negative Rayleigh on Pi_N(S)), as it does for F_theta (S91)?
HIGH-PRECISION NUMERICAL (mpmath dps 90 forms; float64 prolates/angles). Not interval. RH OPEN."""
import os, sys, json, numpy as np, mpmath as mp
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "../s91"); sys.path.insert(0, "../s89")
import h90a as HA, forms as FM
from band_split import delta_ht
from prolate import prolates, ps_eval, chi
from multiprocessing import Pool
N = 24
def span(L):
    lam = np.exp(L); mu = lam**2; s, w = HA.quad_nodes(L, lam); B = HA.basis(L, s)
    M = int(2 * mu) + 12; c, co = prolates(mu, M, K=160)
    chis = [chi(mu, co[n]) for n in range(M)]; nu = max(n for n in range(M) if abs(chis[n]) >= 0.9)
    ps0 = [ps_eval(co[n], 0.0) for n in range(M)]
    def phi(n):
        n0 = n % 2
        return lambda x: ps_eval(co[n], x / lam) * ps0[n0] - ps_eval(co[n0], x / lam) * ps0[n]
    out = {}
    for par in "eo":
        ns = [n for n in range(2, nu + 1) if n % 2 == (0 if par == "e" else 1)]
        So, _ = HA.orth_L2(np.array([HA.Emap(phi(n), s, lam) for n in ns]), w)
        out[par] = (So, (So * w) @ B.T)
    return out, s, w, B
def restricted(Q, Cs):
    CW = mp.matrix(Cs.T.tolist()); Qs = CW.T * Q * CW; Gs = CW.T * CW; Lc = mp.cholesky(Gs) ** -1
    return sorted(mp.eigsy(Lc * Qs * Lc.T, eigvals_only=True))
def task(L):
    mp.mp.dps = 90; sp, s, w, B = span(L); d = FM.load("ZETA", 0, L, N); res = {"L": L}
    for par in "eo":
        Q0 = FM.Qform(d, par); So, Cs = sp[par]; r = {}
        lmax = max(mp.eigsy(Q0, eigvals_only=True))
        r["zeta_on_S_min_rel"] = mp.nstr(restricted(Q0, Cs)[0] / lmax, 5)
        for t in ("-1", "-2", "-4"):
            D, Mz = delta_ht(t, mp.mpf(L), N, par); Q = Q0 + D
            E, V = mp.eigsy(Q); i0 = min(range(N), key=lambda i: E[i])
            g = np.array([float(V[k, i0]) for k in range(N)]); g /= np.linalg.norm(g)
            ev = restricted(Q, Cs)
            r["t" + t] = {"Mz": Mz, "q_min_full_rel": mp.nstr(E[i0] / lmax, 5), "n_neg_full": sum(1 for x in E if x < 0),
                          "min_on_S_rel": mp.nstr(ev[0] / lmax, 5), "n_neg_on_S": sum(1 for x in ev if x < 0), "dimS": len(ev),
                          "sin(ground,S)": float(HA.sin_angles((g @ B)[None, :], So, w)[0])}
        res[par] = r; print(L, par, r, flush=True)
    return res
if __name__ == "__main__":
    Ls = [float(x) for x in sys.argv[1:]] or [0.6, 0.8, 1.0, 1.2]
    with Pool(len(Ls)) as p: out = p.map(task, Ls)
    json.dump(out, open("ht_espan.json", "w"), indent=1)
