"""Grok S94 follow-up: (i) mu_max(C, Q_SS) for zeta vs DH at L 1.3, 1.4 (N 24; DH forms built fresh via s91 forms.load -> s87 hp_forms);
(ii) N = 32 convergence of zeta's mu_max and v^T C v / v^T Q_SS v at L 0.6, 1.0, 1.2. HIGH-PRECISION NUMERICAL (dps 90). RH OPEN."""
import os, sys, json, mpmath as mp
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "../s93"); sys.path.insert(0, "../s92"); sys.path.insert(0, "../s91"); sys.path.insert(0, "../s89"); sys.path.insert(0, "../s90")
import forms as FM, h90a as HA, ht_espan as HE, coupling as CP
from sperp import orthobasis, complement
from multiprocessing import Pool
def task(a):
    L, N, names = a; mp.mp.dps = 90; HA.N = N; HE.N = N; CP.N = N
    sp, s, w, B = HE.span(L); res = {"L": L, "N": N}
    ds = {"ZETA": FM.load("ZETA", 0, L, N)} | ({"DH": FM.load("DH", 0, L, N)} if "DH" in names else {})
    for par in "eo":
        So, Cs = sp[par]; Sb, Sp = complement(Cs); BS = orthobasis(Sb); U = orthobasis(Sp); sg = 1 if par == "e" else -1
        lmax = max(mp.eigsy(FM.Qform(ds["ZETA"], par), eigvals_only=True)); r = {}
        for k, d in ds.items():
            K = mp.zeros(N, N)
            for v in d["p" + par]: K += sg * v * v.T
            o, _ = CP.analyse({"a": d["A" + par], "p": -d["P" + par], "k": K}, None, BS, U, lmax)
            r[k] = {x: o.get(x) for x in ("dimS", "QSS_min_rel", "Qperp_min_rel", "Sigma_min_rel", "Sigma_n_neg", "Sigma_chol_PD", "Sigma_ldl_n_neg", "mu_max(C,QSS)")} | {"vCv/vQv": o.get("along_vmin", {}).get("vCv/vQv")}
            print(L, N, par, k, json.dumps(r[k]), flush=True)
        res[par] = r
    return res
if __name__ == "__main__":
    args = [(1.3, 24, ("ZETA", "DH")), (1.4, 24, ("ZETA", "DH")), (0.6, 32, ("ZETA",)), (1.0, 32, ("ZETA",)), (1.2, 32, ("ZETA",))]
    with Pool(5) as p: out = p.map(task, args)
    json.dump(out, open("dh_extend.json", "w"), indent=1, default=str)
