"""S21 (NUMERICAL): zeta first-prime-window minimizer vs zeta zeros.
For h in a grid, Ritz-minimize Q_zeta (S17 builder, Legendre*(1-u^2)^4 basis, real g) per parity, then:
 (1) zero side 2*sum_{gamma<200} |ghat(gamma)|^2 / ||g||^2 vs Ritz lambda (normalization check),
 (2) real zeros of ghat_min on (0,60) vs zeta ordinates gamma_k (Connes-Consani-Moscovici-type test),
 (3) share of the zero-side sum carried by each of the first zeros."""
import sys, os, json; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s17"))
from s17_build import *
import mpmath as mp, numpy as np
mp.mp.dps = 60
M = int(sys.argv[1]); HS = sys.argv[2:]
Z = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s17", "zeros_T200.json")))["zeta"]
gam = np.array(Z)
def tomp(X, I): return mp.matrix([[mp.mpf(X[i, j].mid().str(70, radius=False)) for j in I] for i in I])
out = []
for hs in HS:
    a, b = (hs.split("/")+["1"])[:2]
    h = arb(a)/arb(b) if "/" in hs else (arb(3).log()/2 if hs == "hstar" else arb(hs))
    hf = float(h.mid().str(30, radius=False))
    B = build(M, h); F = forms(B, h); Q = F["zeta"]; N = B["N"]
    xq, wq = np.polynomial.legendre.leggauss(600); x = xq*hf; w = wq*hf; u = xq
    basis = np.array([np.polynomial.legendre.legval(u, [0]*j+[1])*(1-u**2)**4 for j in range(M)])
    for sector, I in (("even", list(range(0, M, 2))), ("odd", list(range(1, M, 2)))):
        q, n = tomp(Q, I), tomp(N, I); L = mp.cholesky(n); Li = mp.inverse(L)
        ev, V = mp.eigsy(Li*q*Li.T); k = min(range(len(ev)), key=lambda i: ev[i]); v = Li.T*V[:, k]
        c = np.zeros(M); 
        for t, j in enumerate(I): c[j] = float(v[t])
        g = c @ basis; nrm = float(np.sum(w*g*g))
        ghat = lambda t: (g*w) @ np.exp(1j*np.outer(x, np.atleast_1d(t)))
        gz = np.abs(ghat(gam))**2; zs = 2*gz.sum()/nrm
        # real zeros of ghat (real for even g, i*real for odd g)
        ts = np.linspace(0.05, 60, 60000); vals = ghat(ts); f = vals.real if sector == "even" else vals.imag
        roots = [float(ts[i]-f[i]*(ts[i+1]-ts[i])/(f[i+1]-f[i])) for i in range(len(ts)-1) if f[i]*f[i+1] < 0]
        match = [(r, float(gam[np.argmin(abs(gam-r))]), float(abs(gam-r).min())) for r in roots]
        share = (2*gz[:6]/nrm/zs).tolist()
        row = {"h": hs, "M": M, "sector": sector, "ritz": mp.nstr(ev[k], 10), "zero_side_T200": f"{zs:.10e}",
               "ghat_roots_0_60": [round(r, 4) for r in roots], "nearest_gamma_dist": [f"{m[2]:.2e}" for m in match],
               "share_first6": [f"{s:.3f}" for s in share]}
        out.append(row); print(json.dumps(row), flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"out_minimizer_M{M}.json"), "w"), indent=1)
