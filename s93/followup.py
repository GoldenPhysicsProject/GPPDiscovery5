"""Grok S93 follow-up: (1) pole-free / arch-vs-primes on S^perp for zeta (is the pole or Lambda>=0 needed there?);
(2) the complementary Haynsworth split along S (+) S^perp:  Q >= 0 <=> Q_perp > 0 and Sigma_S = Q_SS - Q_S,perp Q_perp^{-1} Q_perp,S >= 0;
    reports min Sigma_S vs min Q|_S for every form (coupling correction);  (3) N = 32 convergence of zeta on S^perp (L 0.6, 1.0, 1.2).
HIGH-PRECISION NUMERICAL (dps 90). Not interval. RH OPEN."""
import os, sys, json, numpy as np, mpmath as mp
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "../s92"); sys.path.insert(0, "../s91"); sys.path.insert(0, "../s89"); sys.path.insert(0, "../s90")
import h90a as HA, forms as FM, ht_espan as HE
from band_split import delta_ht
from sperp import orthobasis, complement, chol_ok
from multiprocessing import Pool
def emin(M): return min(mp.eigsy(M, eigvals_only=True))
def task(arg):
    L, N = arg; mp.mp.dps = 90; HA.N = N; HE.N = N
    sp, s, w, B = HE.span(L); res = {"L": L, "N": N}
    dz = FM.load("ZETA", 0, L, N)
    extra = {} if N != 24 else {"DH": FM.load("DH", 0, L, N), "FT0.1": FM.load("FT", 0.1, L, N)}
    for par in "eo":
        So, Cs = sp[par]; Sb, Sp = complement(Cs); BS = orthobasis(Sb); U = orthobasis(Sp)
        Qz = FM.Qform(dz, par); lmax = max(mp.eigsy(Qz, eigvals_only=True)); A = dz["A" + par]; P = dz["P" + par]
        r = {"zeta_Sperp_min_rel": mp.nstr(emin(U.T * Qz * U) / lmax, 5), "zeta_full_min_rel": mp.nstr(emin(Qz) / lmax, 5)}
        if N == 24:
            Ap = U.T * A * U; Pp = U.T * P * U; R = Ap - Pp
            r["R=A-P_on_Sperp_min_rel"] = mp.nstr(emin(R) / lmax, 5); r["R_on_Sperp_chol"] = chol_ok(R)
            Ea, Va = mp.eigsy(Ap); Ais = Va * mp.diag([1 / mp.sqrt(x) for x in Ea]) * Va.T
            r["rho_perp=lmax(A^-1/2 P A^-1/2)"] = mp.nstr(max(mp.eigsy(Ais * Pp * Ais, eigvals_only=True)), 6)
            r["A_on_Sperp_min_rel"] = mp.nstr(min(Ea) / lmax, 5)
            # Lambda -> |Lambda| sign-free bound: replace P by its absolute value form? report ||P_perp|| / min A_perp
            r["normP_perp/minA_perp"] = mp.nstr(max(abs(x) for x in mp.eigsy(Pp, eigvals_only=True)) / min(Ea), 5)
            forms = {"ZETA": Qz} | {k: FM.Qform(d, par) for k, d in extra.items()}
            for t in ("-2", "-4"):
                D, _ = delta_ht(t, mp.mpf(L), N, par); forms["HT" + t] = Qz + D
            r["split_S_Sperp"] = {}
            for k, Q in forms.items():
                Qss = BS.T * Q * BS; Qsp = BS.T * Q * U; Qpp = U.T * Q * U
                qpp = emin(Qpp)
                o = {"Q|S_min_rel": mp.nstr(emin(Qss) / lmax, 5), "Q_perp_min_rel": mp.nstr(qpp / lmax, 5)}
                if qpp > 0:
                    Sig = Qss - Qsp * mp.inverse(Qpp) * Qsp.T; E = mp.eigsy(Sig, eigvals_only=True)
                    o["Sigma_S_min_rel"] = mp.nstr(min(E) / lmax, 5); o["Sigma_S_n_neg"] = sum(1 for x in E if x < 0)
                    o["Sigma_S_chol_PD"] = chol_ok(Sig) if min(E) > 0 else False
                r["split_S_Sperp"][k] = o
        res[par] = r; print(L, N, par, json.dumps(r), flush=True)
    return res
if __name__ == "__main__":
    args = [(L, 24) for L in (0.6, 0.8, 1.0, 1.2)] + [(L, 32) for L in (0.6, 1.0, 1.2)]
    with Pool(7) as p: out = p.map(task, args)
    json.dump(out, open("followup.json", "w"), indent=1)
