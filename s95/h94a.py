"""Grok S95 -- (1) audit of Turbo wave 229 (S_+ = half-line Schur complement across u = 0) in our basis, and
(2) H94a: is Q_SS = Q(E phi_n, E phi_m) (Connes-Consani prolate vectors, 2106.01715 sec.3) ~ diagonal, with diagonal set by
the Slepian-Pollak concentration defect 1 - |chi_n| ?
(1) Our basis: Legendre P_k(2s/L - 1) on the half-line s in [0, L] (u = log x), extended even (Q_e) / odd (Q_o), same
normalisation.  The half-line functions g_k = (e_k + o_k)/2 (supported on s > 0) and their reflections theta g_k = (e_k - o_k)/2
give A = Q(g, g) = (Q_e + Q_o)/4 and B = Q(g, theta g) = (Q_e - Q_o)/4 (cross-parity terms vanish, CC Lemma 2.5).
S_+ = A - B A^{-1} B.  Exact algebra: S_+ = (A - B) A^{-1} (A + B) = Q_o (Q_e + Q_o)^{-1} Q_e,
A^{-1/2} S_+ A^{-1/2} = I - M^2 with M = A^{-1/2} B A^{-1/2} symmetric, so (A > 0): S_+ >= 0 <=> ||M|| <= 1 <=> Q_e >= 0 and Q_o >= 0,
and n_-(S_+) = n_-(Q_e) + n_-(Q_o).
(2) raw prolate vectors c_n = Pi_N E(phi_n) (float64), Q_SS = C^T Q C, Gram G = C^T C.
HIGH-PRECISION NUMERICAL (mpmath dps 90 forms; float64 prolates). Not interval. RH OPEN."""
import os, sys, json, numpy as np, mpmath as mp
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "../s93"); sys.path.insert(0, "../s92"); sys.path.insert(0, "../s91"); sys.path.insert(0, "../s89")
import h90a as HA, forms as FM
from band_split import delta_ht
from prolate import prolates, ps_eval, chi
from sperp import chol_ok, ldl_inertia
from multiprocessing import Pool
N = 24
nst = lambda x, k=4: mp.nstr(x, k)
def mpm(a): return mp.matrix([[mp.mpf(float(x)) for x in r] for r in a])

def raw_span(L):
    lam = np.exp(L); mu = lam**2; s, w = HA.quad_nodes(L, lam); B = HA.basis(L, s)
    M = int(2 * mu) + 12; c, co = prolates(mu, M, K=160)
    chis = [chi(mu, co[n]) for n in range(M)]; nu = max(n for n in range(M) if abs(chis[n]) >= 0.9)
    ps0 = [ps_eval(co[n], 0.0) for n in range(M)]
    out = {}
    for par in "eo":
        ns = [n for n in range(2, nu + 1) if n % 2 == (0 if par == "e" else 1)]
        rows = []
        for n in ns:
            n0 = n % 2
            f = lambda x, n=n, n0=n0: ps_eval(co[n], x / lam) * ps0[n0] - ps_eval(co[n0], x / lam) * ps0[n]
            g = HA.Emap(f, s, lam); rows.append((g * w) @ B.T)
        out[par] = (ns, np.array(rows), [1 - abs(chis[n]) for n in ns])
    return out

def turbo(Qe, Qo, lmax):
    A = (Qe + Qo) / 4; Bm = (Qe - Qo) / 4
    r = {"A_min_rel": nst(min(mp.eigsy(A, eigvals_only=True)) / lmax), "A_chol_PD": chol_ok(A)}
    Ai = mp.inverse(A); S = A - Bm * Ai * Bm; S = (S + S.T) / 2
    Es = mp.eigsy(S, eigvals_only=True); r["Splus_min_rel"] = nst(min(Es) / lmax, 6); r["Splus_n_neg_eig"] = sum(1 for x in Es if x < 0)
    if min(Es) > 0: r["Splus_chol_PD"] = chol_ok(S)
    else: r["Splus_ldl_n_neg"] = ldl_inertia(S)[0]
    r["identity_err"] = nst(mp.norm(S - Qo * mp.inverse(Qe + Qo) * Qe) / mp.norm(S), 3)
    if r["A_chol_PD"]:
        Lc = mp.cholesky(A); Li = mp.inverse(Lc); Mm = Li * Bm * Li.T
        Em = mp.eigsy((Mm + Mm.T) / 2, eigvals_only=True); r["1-||M||"] = nst(1 - max(abs(x) for x in Em), 6)
    else: r["A_ldl_n_neg"] = ldl_inertia(A)[0]; r["1-||M||"] = "A not PD: S_+ criterion not applicable"
    ee = mp.eigsy(Qe, eigvals_only=True); eo = mp.eigsy(Qo, eigvals_only=True)
    r["qmin_e_rel"] = nst(min(ee) / lmax, 6); r["qmin_o_rel"] = nst(min(eo) / lmax, 6)
    r["n_neg_Qe+Qo"] = sum(1 for x in ee if x < 0) + sum(1 for x in eo if x < 0)
    return r

def h94(Q, C, defect, lmax):
    Cm = mpm(C.T); Qs = Cm.T * Q * Cm; G = Cm.T * Cm; k = Qs.rows
    d = [Qs[i, i] for i in range(k)]; ray = [d[i] / G[i, i] for i in range(k)]
    off = max((abs(Qs[i, j]) / mp.sqrt(abs(d[i] * d[j])) for i in range(k) for j in range(k) if i != j), default=mp.mpf(0))
    # whitened (orthonormal-in-span) spectrum vs sorted defects
    try:
        Lg = mp.inverse(mp.cholesky(G)); W = Lg * Qs * Lg.T; Ew = sorted(mp.eigsy((W + W.T) / 2, eigvals_only=True))
    except ValueError: Ew = []
    lr = [float(mp.log10(abs(x) / lmax)) for x in ray]; ld = [np.log10(max(x, 1e-300)) for x in defect]
    pos = all(x > 0 for x in ray)
    corr = float(np.corrcoef(lr, ld)[0, 1]) if k > 2 else None
    return {"rayleigh_rel": [nst(x / lmax, 3) for x in ray], "defect": [f"{x:.2e}" for x in defect],
            "all_diag_pos": pos, "ray/defect": [nst(ray[i] / lmax / defect[i], 3) for i in range(k)], "max_offdiag_ratio": nst(off, 4), "corr_log(ray)_log(defect)": corr,
            "slope": float(np.polyfit(ld, lr, 1)[0]) if k > 2 else None,
            "whitened_eigs_rel": [nst(x / lmax, 3) for x in Ew]}

def task(L):
    mp.mp.dps = 90; res = {"L": L}; sp = raw_span(L)
    dz = FM.load("ZETA", 0, L, N); dd = FM.load("DH", 0, L, N); df = FM.load("FT", 0.1, L, N)
    Q = {}
    for par in "eo":
        Q[("ZETA", par)] = FM.Qform(dz, par); Q[("DH", par)] = FM.Qform(dd, par); Q[("FT0.1", par)] = FM.Qform(df, par)
        D, _ = delta_ht("-2", mp.mpf(L), N, par); Q[("HT-2", par)] = Q[("ZETA", par)] + D
    lm = {par: max(mp.eigsy(Q[("ZETA", par)], eigvals_only=True)) for par in "eo"}; lmax = max(lm.values())
    res["turbo"] = {F: turbo(Q[(F, "e")], Q[(F, "o")], lmax) for F in ("ZETA", "DH", "FT0.1", "HT-2")}
    res["h94a"] = {}
    for par in "eo":
        ns, C, defect = sp[par]; res["h94a"][par] = {"ns": ns}
        for F in ("ZETA", "DH", "FT0.1", "HT-2"):
            res["h94a"][par][F] = h94(Q[(F, par)], C, defect, lm[par])
    print(json.dumps(res, default=str), flush=True)
    return res

if __name__ == "__main__":
    Ls = [float(x) for x in sys.argv[1:]] or [0.6, 0.8, 1.0, 1.2]
    with Pool(len(Ls)) as p: out = p.map(task, Ls)
    json.dump(out, open("h94a.json", "w"), indent=1, default=str)
