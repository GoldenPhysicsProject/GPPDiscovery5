"""Grok S53: diagnose the most negative eigenvector of the S37 odd Atil (L=27/20,T=700,K=780) dump.
NUMERICAL (inverse power iteration on midpoints, 384 bits). Reports degree profile of the eigenvector
(odd Legendre degrees d = 1+2i) and the spectral location t ~ d/L of its mass."""
import json, sys, math
from flint import arb, arb_mat, ctx
ctx.prec = 384
path = sys.argv[1]; L = float(eval(sys.argv[2])) if len(sys.argv) > 2 else 1.35
A = arb_mat([[arb(x) for x in row] for row in json.load(open(path))]); K = A.nrows()
Ai = A.inv().mid()
v = arb_mat([[arb(1) / (i + 1)] for i in range(K)])
for it in range(80):
    u = Ai * v; nu = sum((u[i, 0] ** 2 for i in range(K)), arb(0)).sqrt()
    v = arb_mat([[(u[i, 0] / nu).mid()] for i in range(K)])
Av = A * v; lam = sum((Av[i, 0] * v[i, 0] for i in range(K)), arb(0))
res = sum(((Av[i, 0] - lam * v[i, 0]) ** 2 for i in range(K)), arb(0)).sqrt()
c2 = [float(v[i, 0]) ** 2 for i in range(K)]
cum = 0; q = {}
for i, c in enumerate(c2):
    cum += c
    for p in (0.1, 0.5, 0.9, 0.99):
        if p not in q and cum >= p: q[p] = 1 + 2 * i
imax = max(range(K), key=lambda i: c2[i])
bands = {}
for i, c in enumerate(c2):
    b = int((1 + 2 * i) / L // 100) * 100; bands[b] = bands.get(b, 0) + c
out = dict(file=path, K=K, lam_rayleigh=float(lam), resid=float(res), peak_degree=1 + 2 * imax, peak_t=(1 + 2 * imax) / L,
           degree_quantiles={str(k): v_ for k, v_ in q.items()}, t_band_mass={str(k): bands[k] for k in sorted(bands)})
print(json.dumps(out, indent=1))
