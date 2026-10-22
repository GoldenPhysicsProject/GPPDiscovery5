"""Grok S93 -- the Weil window form on S^perp, the complement of the Connes-Consani near-radical span S = span E(prolates)
(CC 2106.01715 sec.3; S91 H90a span, nu from |chi| >= 0.9), and the in-band Schur complement restricted to S^perp.
For each form Q (ZETA, DH, FT(0.1), HT(-2), HT(-4) via the s89 zero swap with MZ = 50) and parity, with U an exact (mp, dps 90)
orthonormal basis of the coefficient complement of Pi_N(S) (coefficient metric = L^2(ds), N = 24):
  full      : eig(Q)
  S         : Q compressed to Pi_N(S)                     (S91/S92 quantity)
  Sperp     : Q_perp = U^T Q U
  Sperp_band: Haynsworth split of Q_perp along the arch A_perp = U^T A U: out-of-band block (A_perp >= 0 eigvecs) and the
              in-band Schur complement on the arch-negative (Gamma-band) eigvecs of A_perp.  Q_perp >= 0 <=> out PD and Schur PSD.
  Sperp_low : Q compressed to orth(P_perp V_low), V_low = zeta's eigvecs with eig < 1e-6 lambda_max (fixed reference subspace).
Positivity claims are checked by mp Cholesky (LDL^T) at dps 90 of the compressed matrix; failures by a negative pivot.
Part (c): zeta with one prime removed / sign-flipped (all its powers), Lambda scaled by 1 +- 1% (s90 Pmats atoms).
Interlacing (Cauchy): min over a subspace >= min over everything, so S^perp can never fail at a smaller L than the full form.
HIGH-PRECISION NUMERICAL (mpmath dps 90 forms; float64 prolates/angles define the subspaces). Not interval. RH OPEN."""
import os, sys, json, numpy as np, mpmath as mp
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "../s92"); sys.path.insert(0, "../s91"); sys.path.insert(0, "../s89"); sys.path.insert(0, "../s90")
import h90a as HA, forms as FM
from ht_espan import span
from band_split import delta_ht
from prime_anatomy import Pmats, primepowers
from multiprocessing import Pool
N = 24

def mpm(a): return mp.matrix([[mp.mpf(float(x)) for x in r] for r in a])

def orthobasis(Cn):
    """Cn: float k x N rows spanning a subspace. Returns mp N x m orthonormal basis (Gram-Schmidt twice at dps 90)."""
    M = mpm(np.asarray(Cn).T); cols = []
    for j in range(M.cols):
        v = M[:, j]
        for _ in range(2):
            for c in cols: v = v - (c.T * v)[0] * c
        nv = mp.norm(v)
        if nv > mp.mpf("1e-10"): cols.append(v / nv)
    B = mp.matrix(M.rows, len(cols))
    for j, c in enumerate(cols):
        for r in range(M.rows): B[r, j] = c[r]
    return B

def complement(Cs):
    _, sv, Vt = np.linalg.svd(Cs); k = int(np.sum(sv > 1e-12 * sv[0]))
    return Vt[:k], Vt[k:]

def chol_ok(M, shift=0):
    try:
        mp.cholesky(M - shift * mp.eye(M.rows)); return True
    except (ValueError, ZeroDivisionError):
        return False

def ldl_inertia(M):
    """LDL^T without pivoting at dps 90 -> number of negative pivots (equals n_- if no zero pivot occurs)."""
    n = M.rows; A = M.copy(); neg = 0; dmin = None
    for k in range(n):
        d = A[k, k]
        if d == 0: return None, None
        neg += d < 0; dmin = d if dmin is None or d < dmin else dmin
        for i in range(k + 1, n):
            f = A[i, k] / d
            for j in range(k + 1, i + 1): A[i, j] -= f * A[j, k]; A[j, i] = A[i, j]
    return neg, dmin

def stats(M, lmax, chk=True):
    E = sorted(mp.eigsy(M, eigvals_only=True)); nneg = sum(1 for x in E if x < 0)
    r = {"dim": M.rows, "min_rel": mp.nstr(E[0] / lmax, 5), "n_neg": nneg}
    if chk:
        if E[0] > 0: r["chol_PD"] = chol_ok(M); r["chol_half_margin"] = chol_ok(M, E[0] / 2)
        else:
            ni, _ = ldl_inertia(M); r["ldl_n_neg"] = ni
    return r, E[0]

def analyse(Q, A, B_S, U, W, lmax):
    out = {}
    out["full"], _ = stats(Q, lmax, chk=False)
    out["S"], _ = stats(B_S.T * Q * B_S, lmax, chk=False)
    Qp = U.T * Q * U; out["Sperp"], qp = stats(Qp, lmax)
    Ap = U.T * A * U; EA, VA = mp.eigsy(Ap)
    gi = [i for i in range(Ap.rows) if EA[i] < 0]; pi_ = [i for i in range(Ap.rows) if EA[i] >= 0]
    VG = mp.matrix(Ap.rows, len(gi)); VP = mp.matrix(Ap.rows, len(pi_))
    for c, i in enumerate(gi):
        for r in range(Ap.rows): VG[r, c] = VA[r, i]
    for c, i in enumerate(pi_):
        for r in range(Ap.rows): VP[r, c] = VA[r, i]
    Qpp = VP.T * Qp * VP; o, _ = stats(Qpp, lmax, chk=False); out["Sperp_out"] = o
    if gi:
        Sg = VG.T * Qp * VG - (VG.T * Qp * VP) * mp.inverse(Qpp) * (VP.T * Qp * VG)
        out["Sperp_in_schur"], _ = stats(Sg, lmax, chk=False)
    out["Sperp_nG"] = len(gi)
    out["Sperp_low"], _ = stats(W.T * Q * W, lmax, chk=False) if W.cols else ({"dim": 0}, None)
    return out

def task(L):
    mp.mp.dps = 90; sp, s, w, Bq = span(L); res = {"L": L}
    dz = FM.load("ZETA", 0, L, N); dd = FM.load("DH", 0, L, N); df = FM.load("FT", 0.1, L, N)
    X = mp.exp(2 * mp.mpf(L)); pp = primepowers(X); primes = sorted({p for _, p, _ in pp})
    Pn = {n: Pmats([(n, 1)], mp.mpf(L), N) for n, _, _ in pp}
    for pi_, par in enumerate("eo"):
        So, Cs = sp[par]; Sb, Sp = complement(Cs)
        B_S = orthobasis(Sb); U = orthobasis(Sp)
        Qz = FM.Qform(dz, par); A = dz["A" + par]
        Ez, Vz = mp.eigsy(Qz); lmax = max(Ez)
        low = [i for i in range(N) if Ez[i] < mp.mpf("1e-6") * lmax]
        Vl = np.array([[float(Vz[r, i]) for i in low] for r in range(N)])
        Pp = Sp.T @ Sp; Wl = (Pp @ Vl).T if len(low) else np.zeros((0, N))
        W = orthobasis(Wl) if len(low) else mp.matrix(N, 0)
        # how much of zeta's low eigenspace lies in S^perp (coefficient metric)
        r = {"dim_S_N": B_S.cols, "dim_Sperp": U.cols, "n_low_zeta(<1e-6)": len(low), "dim_Sperp_low": W.cols,
             "low_norm_in_Sperp": [float(np.linalg.norm(Sp @ Vl[:, j])) for j in range(len(low))]}
        Pz = sum((m * Pn[n][pi_] for n, _, m in pp), mp.zeros(N, N))
        EAf = mp.eigsy(A, eigvals_only=True); EAp = mp.eigsy(U.T * A * U, eigvals_only=True)
        r["arch_full_n_neg"] = sum(1 for x in EAf if x < 0); r["arch_on_Sperp_n_neg"] = sum(1 for x in EAp if x < 0)
        r["arch_on_Sperp_min_rel"] = mp.nstr(min(EAp) / lmax, 5); r["arch_full_min_rel"] = mp.nstr(min(EAf) / lmax, 5)
        print(L, par, "dimS", B_S.cols, "dimSperp", U.cols, "nlow", len(low), "arch n- full", r["arch_full_n_neg"], "on Sperp", r["arch_on_Sperp_n_neg"], r["arch_on_Sperp_min_rel"], flush=True)
        r["check_P_rebuild"] = mp.nstr(mp.norm(Pz - dz["P" + par]) / mp.norm(dz["P" + par]), 3)
        forms = {"ZETA": (Qz, A), "DH": (FM.Qform(dd, par), dd["A" + par]), "FT0.1": (FM.Qform(df, par), df["A" + par])}
        for t in ("-2", "-4"):
            D, Mz = delta_ht(t, mp.mpf(L), N, par); forms["HT" + t] = (Qz + D, A)
        sg = 1 if par == "e" else -1
        def Qfrom(P):
            Q = A - P
            for v in dz["p" + par]: Q += sg * v * v.T
            return Q
        for p in primes:
            Pp_ = sum((m * Pn[n][pi_] for n, q, m in pp if q == p), mp.zeros(N, N))
            forms[f"Z_remove_{p}"] = (Qfrom(Pz - Pp_), A); forms[f"Z_flip_{p}"] = (Qfrom(Pz - 2 * Pp_), A)
        for e in ("0.01", "-0.01"): forms[f"Z_scale_{e}"] = (Qfrom((1 + mp.mpf(e)) * Pz), A)
        print(L, par, "P rebuild rel err", r["check_P_rebuild"], flush=True)
        r["forms"] = {}
        for k, (Q, Af) in forms.items():
            r["forms"][k] = analyse(Q, Af, B_S, U, W, lmax)
            f = r["forms"][k]
            print(L, par, k, "full", f["full"]["min_rel"], f["full"]["n_neg"], "| S", f["S"]["min_rel"], f["S"]["n_neg"],
                  "| Sperp", f["Sperp"]["min_rel"], f["Sperp"]["n_neg"], {x: f["Sperp"][x] for x in f["Sperp"] if x.startswith(("chol", "ldl"))},
                  "| out", f["Sperp_out"]["min_rel"], "| inSchur", f.get("Sperp_in_schur", {}).get("min_rel"), "nG", f["Sperp_nG"],
                  "| low", f["Sperp_low"].get("min_rel"), f["Sperp_low"].get("n_neg"), flush=True)
        res[par] = r
    return res

if __name__ == "__main__":
    Ls = [float(x) for x in sys.argv[1:]] or [0.6, 0.8, 1.0, 1.2]
    with Pool(len(Ls)) as p: out = p.map(task, Ls)
    json.dump(out, open("sperp.json", "w"), indent=1, default=str)
