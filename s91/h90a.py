"""Grok S91 -- H90a test: is the near-radical of the Weil window form = E(prolates) (Connes-Consani 2106.01715 sec.3)?
Window: test functions on x in [lambda^-1, lambda], lambda = e^L (our L; S87 builder; primes n <= e^{2L}); CC's mu = lambda^2 = e^{2L}, CC's L_CC = 2L.
Prolates: psi_n(x) = PS_{2n,0}(2 pi mu, x/lambda) on [-lambda, lambda] (time limit lambda, band limit lambda, self-dual Fourier e^{-2 pi i x y}),
  Fourier sign chi_n ~ (-1)^n for n <= nu(mu) ~ 2 mu - 1 (validated against CC's chi(5.5,.), chi(7.5,.) tables to 10 digits).
CC vectors: phi_n = psi_n psi_0(0) - psi_0 psi_n(0) (n even -> even sector), phi_n = psi_n psi_1(0) - psi_1 psi_n(0) (n odd -> odd sector), 2 <= n <= nu;
  g_n(s) = E(phi_n)(e^s) = e^{s/2} sum_{r>=1} phi_n(r e^s) on s = log x in [0, L] (CC eq. 3.4; parity assumed via Poisson).
Our sector basis: b_j(s) = sqrt((2j+1)/L) P_j(2s/L - 1), j < N on [0, L] (parity extension to [-L, L]); coefficient metric = L^2(ds).
All angles are computed in L^2[0, L] with exact quadrature on the pieces between the jumps s = L - log r, so the E-vectors are NOT truncated to N.
Targets (mpmath dps 90): near-null y = R^{-1}p (even, R = A - P), y = R^{-1}q (odd); low eigenspace of Q (eig < 1e-8 lambda_max).
HIGH-PRECISION NUMERICAL for targets, float64 for prolates/angles (floor ~1e-13). Not interval. RH OPEN."""
import os, sys, json, numpy as np, mpmath as mp
os.environ.setdefault("OMP_NUM_THREADS", "1")
from numpy.polynomial import legendre as Lg
from multiprocessing import Pool
from prolate import prolates, ps_eval, chi
import forms as FM
N = 24; NQ = 64; rng = np.random.default_rng(91)

def quad_nodes(L, lam):
    br = sorted({0.0, L} | {L - np.log(r) for r in range(2, int(np.floor(lam)) + 1) if 0 < L - np.log(r) < L})
    x, w = Lg.leggauss(NQ); S, W = [], []
    for a, b in zip(br[:-1], br[1:]):
        S.append((b - a) / 2 * x + (a + b) / 2); W.append((b - a) / 2 * w)
    return np.concatenate(S), np.concatenate(W)

def basis(L, s):
    t = 2 * s / L - 1
    return np.array([np.sqrt((2 * j + 1) / L) * Lg.legval(t, np.eye(N)[j]) for j in range(N)])   # N x nodes

def Emap(fun, s, lam):   # e^{s/2} sum_r f(r e^s), f supported in [-lam, lam]
    out = np.zeros_like(s); x = np.exp(s)
    for r in range(1, int(np.floor(lam)) + 1):
        xr = r * x; m = xr <= lam
        out[m] += fun(xr[m])
    return np.exp(s / 2) * out

def orth_L2(Gfun, w, tol=1e-13):   # Gfun: K x nodes; returns orthonormal function samples (k' x nodes) in L2(w)
    A = Gfun * np.sqrt(w); U, sv, Vt = np.linalg.svd(A.T, full_matrices=False)
    keep = sv > tol * sv[0]
    return (U[:, keep] / np.sqrt(w)[:, None]).T, sv

def sin_angles(Ufun, Sorth, w):   # sin of principal angles of span(Ufun rows, orthonormal) to span(Sorth rows)
    P = (Ufun * w) @ Sorth.T                      # <u_i, s_k>
    R = Ufun - P @ Sorth                          # residuals
    M = R * np.sqrt(w)
    return np.sort(np.linalg.svd(M, compute_uv=False))

def mp_to_np(v): return np.array([float(x) for x in v])

def targets(d, par):
    A = d["A" + par]; P = d["P" + par]; Q = FM.Qform(d, par); R = A - P
    pv = d["p" + par][0] if d["p" + par] else None
    y = mp.lu_solve(R, pv) if pv is not None else None
    E, V = mp.eigsy(Q); idx = sorted(range(N), key=lambda i: E[i]); E = [E[i] for i in idx]
    lmax = E[-1]; low = [k for k in range(N) if E[k] < mp.mpf("1e-8") * lmax]
    Vn = np.array([[float(V[r, idx[k]]) for k in range(N)] for r in range(N)])
    out = {"eigs_low": [mp.nstr(E[k], 6) for k in range(min(8, N))], "lmax": mp.nstr(lmax, 6), "n_low": len(low), "Q": Q, "V": Vn}
    if y is not None:
        yn = mp_to_np(y); out["y"] = yn / np.linalg.norm(yn)
    out["ground"] = Vn[:, 0]
    return out

def task(L):
    mp.mp.dps = 90
    lam = np.exp(L); mu = lam**2
    s, w = quad_nodes(L, lam); B = basis(L, s)
    M = int(2 * mu) + 12
    c, co = prolates(mu, M, K=160)
    chis = [chi(mu, co[n]) for n in range(M)]
    nu = max(n for n in range(M) if abs(chis[n]) >= 0.9)
    nu_strict = max(n for n in range(M) if abs(chis[n]) >= 1 - 1e-6)
    ps0 = [ps_eval(co[n], 0.0) for n in range(M)]
    def phi(n):
        n0 = n % 2
        return lambda x: ps_eval(co[n], x / lam) * ps0[n0] - ps_eval(co[n0], x / lam) * ps0[n]
    G = {n: Emap(phi(n), s, lam) for n in range(2, M)}
    res = {"L": L, "mu=lambda^2": mu, "L_CC=2L": 2 * L, "c=2pi mu": c, "nu(|chi|>=0.9)": nu, "nu_strict(1-|chi|<1e-6)": nu_strict,
           "chi": [round(x, 12) for x in chis], "forms": {}}
    forms = [("ZETA", 0), ("DH", 0), ("FT", 0.1)]
    for F, th in forms:
        d = FM.load(F, th, L, N); rF = {}
        for par in "eo":
            par_n = [n for n in range(2, nu + 1) if n % 2 == (0 if par == "e" else 1)]
            par_ns = [n for n in par_n if n <= nu_strict]
            T = targets(d, par); r = {"eigs_low": T["eigs_low"], "lmax": T["lmax"], "n_low(<1e-8 rel)": T["n_low"], "dim_S(nu)": len(par_n), "dim_S(strict)": len(par_ns)}
            Sorth, sv = orth_L2(np.array([G[n] for n in par_n]), w)
            r["S_rank"] = Sorth.shape[0]
            # fraction of each S direction inside our N=24 polynomial space
            Cs = (Sorth * w) @ B.T                                  # K x N coefficients
            r["S_in_polyN_min_norm2"] = float(np.min(np.sum(Cs**2, 1)))
            # Q restricted to Pi_N(S): generalized Rayleigh (mp)
            CW = mp.matrix(Cs.T.tolist()); Qs = CW.T * T["Q"] * CW; Gs = CW.T * CW
            ev = mp.eigsy(mp.cholesky(Gs)**-1 * Qs * (mp.cholesky(Gs)**-1).T, eigvals_only=True)
            r["Q_on_PiN_S_rel"] = [mp.nstr(x / mp.mpf(T["lmax"]), 4) for x in sorted(ev)]
            fun = lambda v: v @ B
            if "y" in T:
                r["sin(y,S)"] = float(sin_angles(fun(T["y"])[None, :], Sorth, w)[0])
                if par_ns:
                    So2, _ = orth_L2(np.array([G[n] for n in par_ns]), w); r["sin(y,S_strict)"] = float(sin_angles(fun(T["y"])[None, :], So2, w)[0])
            r["sin(ground,S)"] = float(sin_angles(fun(T["ground"])[None, :], Sorth, w)[0])
            nl = T["n_low"]
            if nl:
                r["principal_sin(low,S)"] = [float(x) for x in sin_angles((T["V"][:, :nl].T) @ B, Sorth, w)]
            # matched: first k CC vectors vs k lowest eigvecs
            r["matched"] = []
            for k in range(1, len(par_n) + 1):
                Sk, _ = orth_L2(np.array([G[n] for n in par_n[:k]]), w)
                r["matched"].append({"k": k, "max_sin": float(sin_angles((T["V"][:, :k].T) @ B, Sk, w)[-1]), "eig_k": T["eigs_low"][k - 1] if k <= 8 else None})
            if F == "ZETA":
                K = len(par_n); ctrl = {}
                yv = fun(T["y"])[None, :]
                rs = []
                for _ in range(20):
                    Rr = rng.standard_normal((K, N)) @ B; Ro, _ = orth_L2(Rr, w); rs.append(float(sin_angles(yv, Ro, w)[0]))
                ctrl["gauss_coef_subspace_sin(y) median/min"] = [float(np.median(rs)), float(min(rs))]
                wrong = [n for n in range(nu + 1, M) if n % 2 == (0 if par == "e" else 1)][:K]
                if len(wrong) == K:
                    Wo, _ = orth_L2(np.array([G[n] for n in wrong]), w); ctrl["E(non-concentrated prolates n>nu)_sin(y)"] = float(sin_angles(yv, Wo, w)[0])
                rs = []
                for _ in range(10):
                    Gr = []
                    for k in range(K):
                        a2 = rng.standard_normal(K + 2); cf = np.zeros(2 * K + 6); cf[::2][:K + 2] = a2
                        v0 = Lg.legval(0.0, cf); cf[0] -= v0
                        Gr.append(Emap(lambda x, cf=cf: Lg.legval(x / lam, cf), s, lam))
                    Ro, _ = orth_L2(np.array(Gr), w); rs.append(float(sin_angles(yv, Ro, w)[0]))
                ctrl["E(random even poly f, f(0)=0)_sin(y) median/min"] = [float(np.median(rs)), float(min(rs))]
                r["controls"] = ctrl
            rF[par] = r
            print(L, F, par, {k: v for k, v in r.items() if k not in ("matched", "Q_on_PiN_S_rel")}, "| matched", [(m["k"], "%.2e" % m["max_sin"]) for m in r["matched"]],
                  "| Q|S rel", r["Q_on_PiN_S_rel"][:3], flush=True)
        res["forms"][F + ("" if not th else str(th))] = rF
    return res

if __name__ == "__main__":
    Ls = [float(x) for x in sys.argv[1:]] or [0.6, 0.8, 1.0, 1.2]
    with Pool(len(Ls)) as pool: out = pool.map(task, Ls)
    json.dump(out, open("h90a.json", "w"), indent=1)
