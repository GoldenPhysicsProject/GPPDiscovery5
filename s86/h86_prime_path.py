"""Grok S86 H86a: prime-path (Neumann) expansion of the pole-sewing Schur scalar. NUMERICAL float64. RH OPEN.
Window form on [-L,L], parity split (as s82os/os_hankel.py, S14 normalisation), with kernel multiplier cosh(theta r):
  FORM ZETA theta=0 : zeta.   FORM FT theta>0 : F_theta (zeros rho +- theta; W_F = cosh(theta r) W_zeta, two poles 1/2 +- theta).
  FORM DH : Davenport-Heilbronn (no pole).
Decompose parity form Q = A - P + Pole, A = archimedean+delta (no primes, no pole), P = prime-power part (sum_n Lambda(n)/sqrt n * shift_n),
Pole = +(1/2) sum_alpha p_alpha p_alpha^T (even) / -(1/2) sum q q^T (odd)   [zeta: = p p^T, alpha=1/2].
Pole-free R = A - P.  Neumann: R^{-1} = sum_k A^{-1}(P A^{-1})^k  converges iff rho(A^{-1}P) < 1.
Schur scalar (even, one pole):  s = 1 + p^T R^{-1} p = 1 + sum_k t_k,  t_k = p^T A^{-1}(P A^{-1})^k p  ('k-prime paths').
Two poles: Schur matrix S = I + (1/2)Pm^T R^{-1} Pm (2x2); Weil_even <=> S <= 0 on ... (given n_-(R)=2) -- reported via eigenvalues.
"""
import sys, json, math
import numpy as np
from numpy.polynomial import legendre as leg
sys.path.insert(0, "../s82os")
from os_hankel import mang, dh_weights

def build(FORM, L, N, theta=0.0, M=96, S=160):
    xg, wg = leg.leggauss(M)
    def basis(x):
        x = np.asarray(x, float); t = 2 * x / L - 1
        out = np.zeros((N,) + x.shape); inside = (x >= 0) & (x <= L)
        for k in range(N):
            c = np.zeros(k + 1); c[k] = 1
            out[k] = np.where(inside, math.sqrt((2 * k + 1) / L) * leg.legval(t, c), 0.0)
        return out
    def conv_plus(s):
        a, b = max(0.0, s - L), min(L, s)
        if b <= a: return np.zeros((N, N))
        x = (b - a) / 2 * xg + (a + b) / 2; w = (b - a) / 2 * wg
        return (basis(x) * w) @ basis(s - x).T
    def conv_minus(r):
        a, b = r, L
        if b <= a: return np.zeros((N, N))
        x = (b - a) / 2 * xg + (a + b) / 2; w = (b - a) / 2 * wg
        Cm = (basis(x) * w) @ basis(x - r).T
        return 0.5 * (Cm + Cm.T)
    wt = lambda r: math.cosh(theta * r)
    k1 = lambda z: math.exp(-z) / math.sinh(z)
    k2 = lambda z: 2 * math.exp(-z / 4) * math.sinh(3 * z / 4) / math.sinh(z)
    sech = lambda z: 1 / math.cosh(z / 2)
    x = L / 2 * xg + L / 2; w = L / 2 * wg; P = basis(x)
    mom = lambda al: (P @ (w * np.exp(al * x)), P @ (w * np.exp(-al * x)))
    poles = []
    if FORM == 'ZETA': poles = [(0.5, 2.0)]                     # kernel 2cosh(r/2)
    if FORM == 'FT': poles = [(0.5 + theta, 1.0), (0.5 - theta, 1.0)]  # cosh(r/2)cosh(th r)*2 = cosh((1/2+th)r)+cosh((1/2-th)r)
    pe, po = [], []
    for al, c in poles:   # kernel c*cosh(al r): Toeplitz c/2(a+a-^T+a-a+^T), Hankel c/2(a+a+^T+a-a-^T); even = c/2 (a++a-)(..)^T
        ap, am = mom(al); pe.append(math.sqrt(c / 2) * (ap + am)); po.append(math.sqrt(c / 2) * (ap - am))
    sg, swg = leg.leggauss(S)
    Barch = np.zeros((N, N))
    for (lo, hi) in ((0.0, L), (L, 2 * L)):
        for sk, wk in zip((hi - lo) / 2 * sg + (lo + hi) / 2, (hi - lo) / 2 * swg):
            ker = (-0.5 * (k1(sk) + k2(sk)) + (0.5 * sech(sk) if FORM == 'DH' else 0.0)) * wt(sk)
            Barch += wk * ker * conv_plus(sk)
    E0 = conv_minus(0.0)
    c0 = -0.5772156649015329 - math.log(math.pi) - math.log(1 - math.exp(-2 * L))
    Tarch = c0 * E0 + (math.log(5) * E0 if FORM == 'DH' else 0)
    for sk, wk in zip(L / 2 * sg + L / 2, L / 2 * swg):
        E = conv_minus(sk)
        Tarch += wk * ((E0 - wt(sk) * E) * k1(sk) - wt(sk) * E * k2(sk) + (wt(sk) * E * sech(sk) if FORM == 'DH' else 0))
    NN = int(math.exp(2 * L)) + 1
    wts = dh_weights(NN) if FORM == 'DH' else {n: mang(n) for n in range(2, NN + 1)}
    Tp = np.zeros((N, N)); Bp = np.zeros((N, N))
    for n in range(2, NN + 1):
        if wts.get(n, 0) == 0: continue
        ln = math.log(n); c = wts[n] / math.sqrt(n) * wt(ln)
        if ln < 2 * L: Bp += c * conv_plus(ln)
        if ln < L: Tp += 2 * c * conv_minus(ln)
    sym = lambda X: 0.5 * (X + X.T)
    Tarch, Barch, Tp, Bp = map(sym, (Tarch, Barch, Tp, Bp))
    return dict(Ae=Tarch + Barch, Ao=Tarch - Barch, Pe=Tp + Bp, Po=Tp - Bp,
                pe=np.array(pe).T.reshape(N, -1), po=np.array(po).T.reshape(N, -1))

def analyse(FORM, L, N=28, theta=0.0, kmax=60):
    d = build(FORM, L, N, theta); out = {"form": FORM, "theta": theta, "L": L, "N": N}
    for par, sgn in (("even", +1), ("odd", -1)):
        A, P, V = d["A" + par[0]], d["P" + par[0]], d["p" + par[0]]
        R = A - P; Q = R + sgn * V @ V.T
        eA, eR, eQ = (np.linalg.eigvalsh(X) for X in (A, R, Q))
        Ai = np.linalg.inv(A); K = Ai @ P
        rho = max(abs(np.linalg.eigvals(K)))
        o = {"nneg_A": int((eA < 0).sum()), "nneg_R": int((eR < 0).sum()), "lam_min_Q": eQ[0],
             "rho(A^-1 P)": rho, "A_min": eA[0], "R_min": eR[0]}
        # Cholesky check on Q (positivity statement only if succeeds; float, not a certificate)
        try: np.linalg.cholesky(Q); o["chol_Q"] = "ok(float)"
        except np.linalg.LinAlgError: o["chol_Q"] = "fails"
        if V.shape[1]:
            S = np.eye(V.shape[1]) + sgn * V.T @ np.linalg.solve(R, V)
            o["schur_eigs"] = np.linalg.eigvalsh(0.5 * (S + S.T)).tolist()
            if V.shape[1] == 1:
                v = V[:, 0]; t = []; w_ = Ai @ v
                for k in range(kmax):
                    t.append(float(v @ w_)); w_ = Ai @ (P @ w_)
                o["t_k(first8)"] = t[:8]; o["t0=p'A^-1p"] = t[0]
                o["neumann_partial_s"] = float(1 + sgn * sum(t)) if rho < 1 else "diverges"
                o["sign_pattern_t1..8"] = "".join("+" if x > 0 else "-" for x in t[1:9])
        out[par] = o
    return out

if __name__ == "__main__":
    res = []
    for FORM, th in (("ZETA", 0.0), ("DH", 0.0), ("FT", 0.1), ("FT", 0.25)):
        for L in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7):
            r = analyse(FORM, L, theta=th); res.append(r)
            e, o = r["even"], r["odd"]
            print(FORM, th, L, "| even nA", e["nneg_A"], "nR", e["nneg_R"], "rho %.3g" % e["rho(A^-1 P)"], "lamQ %.3e" % e["lam_min_Q"],
                  "schur", ["%.3e" % s for s in e.get("schur_eigs", [])], e.get("sign_pattern_t1..8", ""),
                  "| odd nA", o["nneg_A"], "nR", o["nneg_R"], "rho %.3g" % o["rho(A^-1 P)"], "lamQ %.3e" % o["lam_min_Q"],
                  "schur", ["%.3e" % s for s in o.get("schur_eigs", [])], o.get("sign_pattern_t1..8", ""), flush=True)
    json.dump(res, open("h86_prime_path_results.json", "w"), indent=1, default=float)
