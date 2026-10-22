"""Grok S78 H77a timing: one cell job of the Mori-route x-space Gram assembly (Cp = V FW, Cp2 = F FW) at full size,
Arb (s71/mori_cert.work_cell, prec 512 + Legendre boost) vs a float64/BLAS version of the same quadrature.
Reports wall time of each and max entrywise |float - Arb mid| (NUMERICAL; the float path has no rigorous bound here).
usage: time_cell.py L N parity cell_index chunk_start chunk_len prec"""
import sys, os, time, json, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's71'))
import mori_cert as mc
from flint import arb, ctx
Ls, N, p, ic, k0, kl, prec = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
ctx.prec = prec; L = mc.parseL(Ls); G = mc.G
G.update(prec=prec, L=L, N=N, p=p); G['nrm'] = [((2 * d + 1) / (2 * L)).sqrt() for d in range(N)]
e2L = (2 * L).exp(); sh = []; n = 2
while arb(n) < e2L:
    lam = mc.mangoldt(n)
    if lam is not None:
        c = arb(n).log() / arb(n).sqrt() if lam == math.log(n) else arb([q for q in range(2, n + 1) if n % q == 0][0]).log() / arb(n).sqrt()
        sh += [(arb(n).log(), c), (-arb(n).log(), c)]
    n += 1
G['shifts'] = sh; cand = []
for (ell, _) in sh:
    if float(ell.mid()) > 0:
        z = ell - L if float((ell - L).mid()) > 0 else L - ell
        if 0 < float(z.mid()) < float(L.mid()) and all(abs(float((z - y).mid())) > 1e-30 for y in cand): cand.append(z)
bpa = [arb(0)] + sorted(cand, key=lambda z: float(z.mid())) + [L]; G['bpa'] = bpa
ks = list(range(k0, k0 + kl)); t0 = time.time()
G['gl'] = {k: arb.legendre_p_root(N, k, weight=True) for k in ks}; tq = time.time() - t0
t0 = time.time(); c1, c2 = mc.work_cell((ic, ks)); ta = time.time() - t0
# float64 path
t0 = time.time()
Lf = float(L.mid()); x0, x1 = float(bpa[ic].mid()), float(bpa[ic + 1].mid()); h = (x1 - x0) / 2; xc = (x0 + x1) / 2
r = np.array([float(G['gl'][k][0].mid()) for k in ks]); w = np.array([float(G['gl'][k][1].mid()) for k in ks])
x = xc + h * r; wk = 2 * h * w; nrm = np.sqrt((2 * np.arange(N) + 1) / (2 * Lf))
def leg(y):  # rows d < N, columns nodes; forward recurrence in double
    P = np.zeros((N, y.size)); P[0] = 1; P[1] = y
    for d in range(1, N - 1): P[d + 1] = ((2 * d + 1) * y * P[d] - d * P[d - 1]) / (d + 1)
    return P * nrm[:, None]
V = leg(x / Lf); F = np.zeros_like(V)
for (a, ca) in sh:
    y = x - float(a.mid()); msk = (y > -Lf) & (y < Lf)
    if msk.any():
        Pm = leg(np.where(msk, y, 0) / Lf); F += float(ca.mid()) * Pm * msk[None, :]
degs = list(range(p, N, 2)); V = V[degs]; F = F[degs]
FW = (F * wk[None, :]).T; Cp = V @ FW; Cp2 = F @ FW; tf = time.time() - t0
dev1 = max(abs(float(arb(c1[i][j]).mid()) - Cp[i, j]) for i in range(len(degs)) for j in range(len(degs)))
dev2 = max(abs(float(arb(c2[i][j]).mid()) - Cp2[i, j]) for i in range(len(degs)) for j in range(len(degs)))
print(json.dumps(dict(L=Ls, N=N, parity=p, cell=ic, nodes=kl, prec=prec, t_roots=round(tq, 2), t_arb=round(ta, 2), t_float=round(tf, 3),
      speedup=round(ta / tf, 1), maxdev_Cp=dev1, maxdev_Cp2=dev2, max_Cp=float(np.abs(Cp).max()))), flush=True)
