"""Grok S94 -- H93a: anatomy of the S-S^perp coupling C = Q_S,perp Q_perp^{-1} Q_perp,S on the Connes-Consani span S
(CC 2106.01715 sec.3; S91-S93 span), so that Sigma_S = Q_SS - C (Haynsworth; Q >= 0 <=> Q_perp > 0 and Sigma_S >= 0).
Forms (s93 conventions, dps 90, N 24): Q = A - P + sg*sum pp^T (+ D for H_-t zero swap, MZ 50).
Off-diagonal block X = B_S^T Q U is LINEAR in the sources:  X = X_a (arch) + X_p (= -B_S^T P U, primes) + X_k (pole) [+ X_d (zero swap)].
With G = Q_perp^{-1} (the form's own), C = sum_{i,j} X_i G X_j^T: diagonal pieces C_ii = X_i G X_i^T are PSD automatically (G > 0),
cross pieces C_ij + C_ji are indefinite in general.  Reported per piece: eigen range, and value along v = unit min eigvec of Sigma_S
(and along v_zeta, zeta's Sigma_S min eigvec) relative to v^T Q_SS v; generalized eigenvalues mu of (C, Q_SS) when Q_SS > 0
(Sigma_S >= 0 <=> mu_max <= 1).  Positivity of Sigma_S: mp Cholesky (also at half margin); failures: LDL^T negative pivots.
Note: G is form-dependent, so 'piece' = linear split of X only (the inverse is not split).
HIGH-PRECISION NUMERICAL (mpmath dps 90; float64 prolates define S). Not interval. RH OPEN."""
import os, sys, json, numpy as np, mpmath as mp
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "../s93"); sys.path.insert(0, "../s92"); sys.path.insert(0, "../s91"); sys.path.insert(0, "../s89"); sys.path.insert(0, "../s90")
import forms as FM
from ht_espan import span
from band_split import delta_ht
from prime_anatomy import Pmats, primepowers
from sperp import orthobasis, complement, chol_ok, ldl_inertia
from multiprocessing import Pool
N = 24
def eig(M): return mp.eigsy(M)
def nst(x, k=4): return mp.nstr(x, k)

def gen_eigs(C, Qss):
    """mu solving C v = mu Qss v (Qss > 0) via Cholesky whitening."""
    Lc = mp.cholesky(Qss); Li = mp.inverse(Lc)
    return sorted(mp.eigsy(Li * C * Li.T, eigvals_only=True))

def analyse(pieces, Gfull, BS, U, lmax, vref=None):
    """pieces: dict name -> N x N matrix whose sum is Q. Returns stats."""
    Q = sum(pieces.values(), mp.zeros(N, N))
    Qss = BS.T * Q * BS; Qpp = U.T * Q * U
    Epp = mp.eigsy(Qpp, eigvals_only=True)
    out = {"Qperp_min_rel": nst(min(Epp) / lmax)}
    if min(Epp) <= 0: out["note"] = "Q_perp not PD; split invalid"; return out, None
    G = mp.inverse(Qpp)
    X = {k: BS.T * M * U for k, M in pieces.items()}
    C = sum((X[k] for k in X), mp.zeros(BS.cols, U.cols)); C = C * G * C.T
    Sig = Qss - C
    Es, Vs = mp.eigsy(Sig); i0 = min(range(Sig.rows), key=lambda i: Es[i]); v = Vs[:, i0]
    Eq = mp.eigsy(Qss, eigvals_only=True); Ec = mp.eigsy(C, eigvals_only=True)
    out.update({"dimS": BS.cols, "QSS_min_rel": nst(min(Eq) / lmax), "QSS_max_rel": nst(max(Eq) / lmax),
                "C_eigs_rel": [nst(x / lmax, 3) for x in sorted(Ec)], "Sigma_min_rel": nst(Es[i0] / lmax),
                "Sigma_n_neg": sum(1 for x in Es if x < 0)})
    if Es[i0] > 0: out["Sigma_chol_PD"] = chol_ok(Sig); out["Sigma_chol_half"] = chol_ok(Sig, Es[i0] / 2)
    else: out["Sigma_ldl_n_neg"] = ldl_inertia(Sig)[0]
    if min(Eq) > 0:
        mu = gen_eigs(C, Qss); out["mu_max(C,QSS)"] = nst(mu[-1], 12); out["1-mu_max"] = nst(1 - mu[-1], 5)
    else: out["QSS_n_neg"] = sum(1 for x in Eq if x < 0)
    def along(w):
        qv = (w.T * Qss * w)[0]; cv = (w.T * C * w)[0]; r = {"vQv_rel": nst(qv / lmax), "vCv/vQv": nst(cv / qv, 10)}
        ks = list(X)
        for a in range(len(ks)):
            for b in range(a, len(ks)):
                M = X[ks[a]] * G * X[ks[b]].T
                if a != b: M = M + M.T
                r[f"{ks[a]}{ks[b]}/vQv"] = nst((w.T * M * w)[0] / qv, 5)
        # pieces of Q_SS itself along w
        for k in ks: r[f"QSS_{k}/vQv"] = nst((w.T * (BS.T * pieces[k] * BS) * w)[0] / qv, 5)
        return r
    out["along_vmin"] = along(v)
    if vref is not None: out["along_vzeta"] = along(vref)
    # sign structure of the cross pieces (full eigen range on S)
    ks = list(X); sg = {}
    for a in range(len(ks)):
        for b in range(a, len(ks)):
            M = X[ks[a]] * G * X[ks[b]].T
            if a != b: M = M + M.T
            E = mp.eigsy(M, eigvals_only=True); sg[ks[a] + ks[b]] = [nst(min(E) / lmax, 3), nst(max(E) / lmax, 3)]
    out["piece_eig_range_rel"] = sg
    chk = sum((X[ks[a]] * G * X[ks[b]].T for a in range(len(ks)) for b in range(len(ks))), mp.zeros(BS.cols, BS.cols))
    out["piece_sum_err"] = nst(mp.norm(chk - C) / (mp.norm(C) + mp.mpf(10) ** -80), 3)
    return out, v

def task(L):
    mp.mp.dps = 90; sp, s, w, Bq = span(L); res = {"L": L}
    dz = FM.load("ZETA", 0, L, N); dd = FM.load("DH", 0, L, N); df = FM.load("FT", 0.1, L, N)
    X_ = mp.exp(2 * mp.mpf(L)); pp = primepowers(X_); primes = sorted({p for _, p, _ in pp})
    Pn = {n: Pmats([(n, 1)], mp.mpf(L), N) for n, _, _ in pp}
    for pi_, par in enumerate("eo"):
        So, Cs = sp[par]; Sb, Sp = complement(Cs); BS = orthobasis(Sb); U = orthobasis(Sp)
        sg = 1 if par == "e" else -1
        def pole(d):
            K = mp.zeros(N, N)
            for v in d["p" + par]: K += sg * v * v.T
            return K
        def base(d, P=None): return {"a": d["A" + par], "p": -(d["P" + par] if P is None else P), "k": pole(d)}
        Qz = FM.Qform(dz, par); lmax = max(mp.eigsy(Qz, eigvals_only=True))
        forms = {"ZETA": base(dz), "DH": base(dd), "FT0.1": base(df)}
        for t in ("-2", "-4"):
            D, _ = delta_ht(t, mp.mpf(L), N, par); forms["HT" + t] = base(dz) | {"d": D}
        Pz = sum((m * Pn[n][pi_] for n, _, m in pp), mp.zeros(N, N))
        for p in primes:
            if p > 7: continue
            Pp_ = sum((m * Pn[n][pi_] for n, q, m in pp if q == p), mp.zeros(N, N))
            forms[f"Z_remove_{p}"] = base(dz, Pz - Pp_); forms[f"Z_flip_{p}"] = base(dz, Pz - 2 * Pp_)
        for e in ("0.01", "-0.01"): forms[f"Z_scale_{e}"] = base(dz, (1 + mp.mpf(e)) * Pz)
        r = {"P_rebuild_err": nst(mp.norm(Pz - dz["P" + par]) / mp.norm(dz["P" + par]), 3), "forms": {}}
        oz, vz = analyse(forms["ZETA"], None, BS, U, lmax); r["forms"]["ZETA"] = oz
        for k, pcs in forms.items():
            if k == "ZETA": continue
            o, _ = analyse(pcs, None, BS, U, lmax, vref=vz); r["forms"][k] = o
        for k, o in r["forms"].items():
            a = o.get("along_vmin", {})
            print(L, par, k, "Sig", o.get("Sigma_min_rel"), o.get("Sigma_n_neg"), "chol", o.get("Sigma_chol_PD"), "ldl", o.get("Sigma_ldl_n_neg"),
                  "| 1-mu", o.get("1-mu_max"), "| vCv/vQv", a.get("vCv/vQv"), "| aa", a.get("aa/vQv"), "pp", a.get("pp/vQv"), "kk", a.get("kk/vQv"),
                  "ap", a.get("ap/vQv"), "ak", a.get("ak/vQv"), "pk", a.get("pk/vQv"), "dd", a.get("dd/vQv"),
                  "| at vz:", o.get("along_vzeta", {}).get("vCv/vQv"), "err", o.get("piece_sum_err"), flush=True)
        res[par] = r
    return res

if __name__ == "__main__":
    Ls = [float(x) for x in sys.argv[1:]] or [0.6, 0.8, 1.0, 1.2]
    with Pool(len(Ls)) as p: out = p.map(task, Ls)
    json.dump(out, open("coupling.json", "w"), indent=1, default=str)
