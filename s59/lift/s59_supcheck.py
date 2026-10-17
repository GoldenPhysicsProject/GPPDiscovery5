"""Grok S59 H59a: independent, heavier 6-D sup ESTIMATE of F = B - m for a torus lift json (K, a_lp or a_q), L = 7/5.
Candidates: products of the top local maxima of the two group-local parts (fine 3-D grids NG^3, backtracking ascent), plus NRAND
uniform points (new seed); top NTOP refined by backtracking gradient ascent on the full 6-D F. Lower estimate of the sup (NUMERICAL).
usage: s59_supcheck.py lift.json [seed]. RH OPEN."""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 's28'))
from s28_torus_lp import setup
d = json.load(open(sys.argv[1])); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 7
K = np.array(d['K'], float); a = np.array([float(x) for x in (d.get('a_lp') or d['a_q'])])
NG = int(os.environ.get('NG', '80')); NR = int(os.environ.get('NRAND', '400000')); NTOP = int(os.environ.get('NTOP', '8000')); CH = 3000
ps, lp, vis = setup(1.4, [2, 3, 5, 7, 11, 13]); rs = np.random.default_rng(seed); t0 = time.time()
cv = np.array([c for (_, _, c) in vis]); iv = np.array([i for (i, _, _) in vis]); kv = np.array([k for (_, k, _) in vis], float)
def F(th, msk=None, vm=None):
    KK = K if msk is None else K[msk]; aa = a if msk is None else a[msk]
    ii, kk, cc = (iv, kv, cv) if vm is None else (iv[vm], kv[vm], cv[vm])
    out = np.empty(len(th))
    for s in range(0, len(th), CH):
        t = th[s:s + CH]; out[s:s + CH] = np.cos(t[:, ii] * kk) @ cc - np.cos(t @ KK.T) @ aa
    return out
def ascend(P, f, fun, iters, h=1e-6):
    lr = np.full(len(P), 0.005)
    for _ in range(iters):
        g = np.zeros_like(P)
        for j in range(6):
            e = np.zeros(6); e[j] = h; g[:, j] = (fun(P + e) - fun(P - e)) / (2 * h)
        nP = P + lr[:, None] * g; fn = fun(nP); ok = fn > f
        P[ok] = nP[ok]; f[ok] = fn[ok]; lr = np.where(ok, lr * 1.3, lr * 0.4)
    return P, f
def gmax(idx, top=150):
    g1 = (np.arange(NG) + 0.5) * 2 * math.pi / NG - math.pi
    G = np.stack(np.meshgrid(g1, g1, g1, indexing='ij'), -1).reshape(-1, 3); th = np.zeros((len(G), 6)); th[:, idx] = G
    msk = np.array([all(k[j] == 0 for j in range(6) if j not in idx) for k in K]); vm = np.isin(iv, idx)
    fun = lambda t: F(t, msk, vm); f = fun(th); o = np.argsort(-f)[:top * 6]; P, f = ascend(th[o].copy(), f[o].copy(), fun, 150)
    o = np.argsort(-f); return P[o][:top], f[o][:top]
PA, fA = gmax([0, 1, 2]); PB, fB = gmax([3, 4, 5])
prod = (PA[:, None, :] + PB[None, :, :]).reshape(-1, 6)
cand = np.vstack([prod, rs.uniform(-math.pi, math.pi, (NR, 6))]); f = F(cand); o = np.argsort(-f)[:NTOP]
P, f = ascend(cand[o].copy(), f[o].copy(), F, 200); o = np.argsort(-f)
print(json.dumps(dict(file=sys.argv[1], seed=seed, NG=NG, NRAND=NR, NTOP=NTOP, groupA_local_max=float(fA[0]), groupB_local_max=float(fB[0]),
      sup_est=float(f[o[0]]), top5=[round(float(x), 6) for x in f[o[:5]]], argmax=P[o[0]].round(5).tolist(), wall_s=round(time.time() - t0))), flush=True)
