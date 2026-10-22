"""Grok S101, H100a: signed prime-sum split on zeta's lambda_2(R_e) vector (Queen step (c) eigenvectors).
For an even window function f on [-L,L] (half-width L), g = f * f~ (autocorrelation on [-2L,2L]).
Prime side P(f) = sum_n w(n) n^{-1/2} 2 g(log n), w = Lambda (zeta) or the DH coefficients (s87 dh_weights).
Split P = P_+ + P_- by sign of g(log n); cancellation ratio rho = |P| / (P_+ + |P_-|).
Controls: zeta lambda_1 (nodeless soft mode), DH's own lambda_2(Q_e) vector, random even unit vectors.
NUMERICAL (float64 quadrature of g; mpmath only for the forms). RH OPEN."""
import sys, json, numpy as np, mpmath as mp
from numpy.polynomial import legendre as Lg
sys.path.insert(0, "s87"); import hp_forms as H
xg, wg = Lg.leggauss(200)
def fvals(c, L, s):
    a = np.abs(s); t = 2 * a / L - 1
    return sum(ci * np.sqrt((2 * j + 1) / L) * Lg.legval(t, [0] * j + [1]) for j, ci in enumerate(c))
def g(c, L, u):
    pts = sorted({-L, L - u, max(-L, -u), 0.0} & set(p for p in [-L, L - u, max(-L, -u), 0.0] if -L <= p <= L - u))
    tot = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        if b - a < 1e-15: continue
        s = (b - a) / 2 * xg + (a + b) / 2
        tot += (b - a) / 2 * np.sum(wg * fvals(c, L, s) * fvals(c, L, s + u))
    return tot
def split(c, L, w):
    c = np.array(c, float); c = c / np.linalg.norm(c) / np.sqrt(2)  # unit on [-L,L]
    terms = []
    for n, wn in w.items():
        if wn == 0 or np.log(n) >= 2 * L: continue
        terms.append(float(wn) / np.sqrt(n) * 2 * g(c, L, np.log(n)))
    t = np.array(terms); Pp = t[t > 0].sum(); Pm = t[t < 0].sum()
    return dict(P=Pp + Pm, Pplus=Pp, Pminus=Pm, rho=abs(Pp + Pm) / (Pp - Pm), nterms=len(t))
def lam_vec(form, L, N, k):
    H.setup(60); H._GN.clear(); d = H.build(form, mp.mpf(L), N)
    M = mp.matrix(H.add(d["Ae"], d["Pe"], 1, -1))
    E, V = mp.eigsy(M); idx = sorted(range(N), key=lambda i: E[i]); i = idx[k]
    return float(E[i]), [float(V[j, i]) for j in range(N)]
out = []; rng = np.random.default_rng(101)
Q = json.load(open("s101/queen_lambda2_Re_eigenvectors.json"))["vectors"]
for rec in Q:
    L = float(rec["L"])
    if L <= 1.4: continue
    c = [float(x) for x in eval(rec["coefficients"])] if isinstance(rec["coefficients"], str) else [float(x) for x in rec["coefficients"]]
    NN = int(np.exp(2 * L)) + 1; wz = {n: float(H.mang(n)) for n in range(2, NN + 1)}; wd = {n: float(v) for n, v in H.dh_weights(NN).items()}
    r = dict(L=L, zeta_lam2=split(c, L, wz))
    l1, v1 = lam_vec('ZETA', L, 24, 0); r["zeta_lam1"] = split(v1, L, wz); r["zeta_lam1_val"] = l1
    l2d, v2d = lam_vec('DH', L, 24, 0); r["DH_lam1Q"] = split(v2d, L, wd); r["DH_lam1Q_val"] = l2d
    l2d2, v2d2 = lam_vec('DH', L, 24, 1); r["DH_lam2Q"] = split(v2d2, L, wd)
    rr = [split(rng.standard_normal(24) / np.arange(1, 25), L, wz)["rho"] for _ in range(20)]
    r["random_rho_median"] = float(np.median(rr))
    out.append(r)
    print(L, {k: (round(v["rho"], 6) if isinstance(v, dict) else v) for k, v in r.items()}, flush=True)
    json.dump(out, open("s101/h100a_hiL.json", "w"), indent=1)
