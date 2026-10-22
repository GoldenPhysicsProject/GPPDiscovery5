"""Grok S96 -- H95a: ||M|| <= 1 as an L -> 2L doubling step.
Basis (S95): orthonormal Legendre on the half-line s in [0, L], even (Q_e) / odd (Q_o) extensions.
A = (Q_e+Q_o)/4 = Q on functions supported in [0, L] (a window of length L, i.e. the full form at half-width L/2 up to span),
B = (Q_e-Q_o)/4, M = A^{-1/2} B A^{-1/2} (Cholesky whitening).  A^{-1/2}(Q_e/2)A^{-1/2} = I + M,  A^{-1/2}(Q_o/2)A^{-1/2} = I - M.
So 1 + mu_min(M) <-> Q_e positivity, 1 - mu_max(M) <-> Q_o positivity (generalised eigenvalues of Q_e/2, Q_o/2 against A).
Per (form, L): lambda_min(A), spectrum of M, extremal eigvec vs Q_e/Q_o near-null vectors, mp Cholesky / LDL^T signs (dps 90).
Doubling: compare A(L) with the full form at L/2 (n_-, lambda_min), and 1-||M(L)|| with q_min(L)/lambda_min(A(L)) and q_min(L/2).
HIGH-PRECISION NUMERICAL (mpmath dps 90), not interval. N=24 values < ~1e-26 rel at L >= 1 may be unconverged (S91). RH OPEN."""
import os, sys, json, mpmath as mp
os.environ.setdefault("OMP_NUM_THREADS", "1")
for d in ("../s93", "../s89", "../s91"): sys.path.insert(0, d)
import forms as FM
from band_split import delta_ht
from sperp import chol_ok, ldl_inertia
from multiprocessing import Pool
N = 24
nst = lambda x, k=5: mp.nstr(x, k)
FORMS = ("ZETA", "DH", "FT0.1", "HT-2", "HT-4")

def sign_cert(Mx):
    if chol_ok(Mx): return {"PD_chol": True, "n_neg": 0}
    n, _ = ldl_inertia(Mx); E = mp.eigsy(Mx, eigvals_only=True)
    return {"PD_chol": False, "n_neg_ldl": n, "n_neg_eig": sum(1 for x in E if x < 0)}

def mineig(Mx):
    E, V = mp.eigsy(Mx); i = min(range(Mx.rows), key=lambda j: E[j]); return E[i], V[:, i]

def unit(v): return v / mp.norm(v)

def analyse(Qe, Qo, lmax):
    A = (Qe + Qo) / 4; Bm = (Qe - Qo) / 4
    qe, ve = mineig(Qe); qo, vo = mineig(Qo); a, va = mineig(A)
    r = {"qmin_e_rel": nst(qe / lmax), "qmin_o_rel": nst(qo / lmax), "Qe": sign_cert(Qe), "Qo": sign_cert(Qo),
         "Amin_rel": nst(a / lmax), "A": sign_cert(A)}
    r["_num"] = {"qe": qe / lmax, "qo": qo / lmax, "a": a / lmax}
    if not r["A"]["PD_chol"]:
        r["M"] = "A not PD: doubling step undefined"; return r
    Lc = mp.cholesky(A); Li = mp.inverse(Lc); M = Li * Bm * Li.T; M = (M + M.T) / 2
    E, V = mp.eigsy(M); idx = sorted(range(M.rows), key=lambda j: E[j])
    mmin, mmax = E[idx[0]], E[idx[-1]]
    r["mu_min"] = nst(mmin, 12); r["mu_max"] = nst(mmax, 12)
    r["1+mu_min (Q_e)"] = nst(1 + mmin); r["1-mu_max (Q_o)"] = nst(1 - mmax)
    r["1-||M||"] = nst(1 - max(abs(mmin), abs(mmax))); r["extremal_side"] = "Q_e (mu~-1)" if abs(mmin) >= abs(mmax) else "Q_o (mu~+1)"
    r["n_|mu|>1"] = sum(1 for x in E if abs(x) > 1)
    r["n_mu<-1"] = sum(1 for x in E if x < -1); r["n_mu>+1"] = sum(1 for x in E if x > 1)
    # extremal eigvectors in coefficient space (orthonormal basis => Euclidean = L^2)
    xe = unit(Li.T * V[:, idx[0]]); xo = unit(Li.T * V[:, idx[-1]])
    r["cos(x_ext_e, v_Qe_min)"] = nst(abs((xe.T * ve)[0]), 8); r["cos(x_ext_o, v_Qo_min)"] = nst(abs((xo.T * vo)[0]), 8)
    r["cos(x_ext_e, v_Amin)"] = nst(abs((xe.T * va)[0]), 4)
    # S_+ near-null vector = extremal vector (exact algebra); check
    S = A - Bm * mp.inverse(A) * Bm; S = (S + S.T) / 2; s, vs = mineig(S)
    xm = xe if abs(mmin) >= abs(mmax) else xo
    r["cos(x_ext, v_Splus_min)"] = nst(abs((xm.T * vs)[0]), 8)
    # Rayleigh identity 1+mu_min = min (x'Q_e x/2)/(x'A x); law check: (q_e/2)/(v_e'A v_e)
    ra = (ve.T * A * ve)[0]
    r["(qe/2)/(ve'A ve) / (1+mu_min)"] = nst((qe / 2 / ra) / (1 + mmin), 8)
    r["(qe/2)/lamminA / (1+mu_min)"] = nst((qe / 2 / a) / (1 + mmin), 5)
    r["ve'A ve / lamminA"] = nst(ra / a, 5)
    r["_num"].update({"onem": 1 - max(abs(mmin), abs(mmax)), "pe": 1 + mmin, "po": 1 - mmax})
    return r

def task(L):
    mp.mp.dps = 90; Lm = mp.mpf(L)
    dz = FM.load("ZETA", 0, L, N); dd = FM.load("DH", 0, L, N); df = FM.load("FT", 0.1, L, N)
    Q = {}
    for par in "eo":
        Q[("ZETA", par)] = FM.Qform(dz, par); Q[("DH", par)] = FM.Qform(dd, par); Q[("FT0.1", par)] = FM.Qform(df, par)
        for t in ("-2", "-4"):
            D, _ = delta_ht(t, Lm, N, par); Q[("HT" + t, par)] = Q[("ZETA", par)] + D
    lmax = max(max(mp.eigsy(Q[("ZETA", p)], eigvals_only=True)) for p in "eo")
    res = {"L": L, "lmax_zeta": nst(lmax), "forms": {F: analyse(Q[(F, "e")], Q[(F, "o")], lmax) for F in FORMS}}
    for F in FORMS: res["forms"][F]["_num"] = {k: str(v) for k, v in res["forms"][F]["_num"].items()}
    res["_lmax"] = str(lmax)
    print(json.dumps({"L": L, "done": True}), flush=True)
    return res

if __name__ == "__main__":
    Ls = [float(x) for x in sys.argv[1:]] or [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
    with Pool(6) as p: out = p.map(task, Ls, chunksize=1)
    json.dump(out, open("h95a.json", "w"), indent=1, default=str)
