"""Grok S59 H59a (test of the S58 H58c conjecture): CROSS-GROUP buffered torus lift at L = 7/5 on the 6-torus of {2,3,5,7,11,13}.
Start = the S58 CERTIFIED separable lift ({2,3,5} D9 + {7,11,13} D8, A_cert <= 5.189). Add cross characters k = (k_A, k_B) with
k_A != 0 != k_B, |k|_1 <= DC, |k.log p| >= 2L + EPS (buffer, EPS = 0.0056 as in S55/S58), and re-optimise ALL coefficients with the
S55 proximal (trust-region) cutting-plane LP. True sup on the 6-torus is ESTIMATED (product of group local maxima + random + previous
active points, then per-point backtracking gradient ascent). NUMERICAL design only (a 6-D sup estimate can be low; rigour would need
B&B). Memory kept < ~1 GB (chunked evaluation). usage: s59_joint_lp.py DC budget_s tag [init.json]. RH OPEN."""
import sys, os, json, math, time, itertools
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 's28'))
from s28_torus_lp import setup
from scipy.optimize import linprog
EPS = float(os.environ.get('EPS', '0.0056')); L = 1.4
DC = int(sys.argv[1]); budget = float(sys.argv[2]); tag = sys.argv[3]; init = sys.argv[4] if len(sys.argv) > 4 else None
NR = int(os.environ.get('NRAND', '60000')); NTOP = int(os.environ.get('NTOP', '3000')); PMAX = int(os.environ.get('PMAX', '5000')); CH = 3000
SEP = '/workspace/GPPDiscovery5/s58/run140/torus_L1.40_buf_D9D8_cert.json'
rs = np.random.default_rng(int(os.environ.get('SEED', '1'))); t0 = time.time()
ps, lp, vis = setup(L, [2, 3, 5, 7, 11, 13]); n = 6
sep = json.load(open(SEP)); K0 = np.array(sep['K'], int); a0 = np.array([float(x) for x in sep['a_q']])
# cross characters
cross = []
for k in itertools.product(range(-DC, DC + 1), repeat=6):
    k = np.array(k)
    if np.abs(k).sum() > DC or not k[:3].any() or not k[3:].any(): continue
    if k[k != 0][0] < 0: continue
    if abs(k @ lp) >= 2 * L + EPS: cross.append(k)
cross = np.array(cross, int).reshape(-1, 6)
K = np.vstack([K0, cross]) if len(cross) else K0; a = np.r_[a0, np.zeros(len(cross))]
if init:
    d = json.load(open(init)); m = {tuple(k): x for k, x in zip(d['K'], d['a_lp'])}
    a = np.array([m.get(tuple(k), 0.0) for k in K.tolist()])
Kf = K.astype(float)
print(f'L={L} EPS={EPS} DC={DC} sep_chars={len(K0)} cross_chars={len(cross)} minfreq={np.abs(K@lp).min():.6f}', flush=True)
cvis = np.array([c for (_, _, c) in vis]); ivis = np.array([i for (i, _, _) in vis]); kvis = np.array([k for (_, k, _) in vis], float)
def Bv(th): return np.cos(th[:, ivis] * kvis) @ cvis
def F(th, a):
    out = np.empty(len(th))
    for s in range(0, len(th), CH):
        t = th[s:s + CH]; out[s:s + CH] = Bv(t) - np.cos(t @ Kf.T) @ a
    return out
def G(th, a):
    g = np.zeros_like(th)
    for s in range(0, len(th), CH):
        t = th[s:s + CH]; S = np.sin(t[:, ivis] * kvis) * (cvis * kvis)
        gg = np.zeros_like(t)
        for j in range(len(ivis)): gg[:, ivis[j]] -= S[:, j]
        gg += (np.sin(t @ Kf.T) * a) @ Kf; g[s:s + CH] = gg
    return g
def refine(th, a, iters=120):
    th = th.copy(); f = F(th, a); lr = np.full(len(th), 0.005)
    for _ in range(iters):
        g = G(th, a); nt = th + lr[:, None] * g; fn = F(nt, a); ok = fn > f
        th[ok] = nt[ok]; f[ok] = fn[ok]; lr = np.where(ok, lr * 1.3, lr * 0.4)
    return th, f
# group local maxima of the separable start (3-D grids), for product candidates
def group_max(idx, NG=56, top=120):
    g1 = (np.arange(NG) + 0.5) * 2 * math.pi / NG - math.pi
    Gd = np.stack(np.meshgrid(g1, g1, g1, indexing='ij'), -1).reshape(-1, 3)
    th = np.zeros((len(Gd), 6)); th[:, idx] = Gd
    # zero out the other group: evaluate group-local part only
    msk = np.array([all(k[j] == 0 for j in range(6) if j not in idx) for k in K0])
    vm = np.isin(ivis, idx)
    def Fg(t, a=a0[msk]):
        return np.cos(t[:, ivis[vm]] * kvis[vm]) @ cvis[vm] - np.cos(t @ K0[msk].astype(float).T) @ a
    f = np.concatenate([Fg(th[s:s + CH]) for s in range(0, len(th), CH)])
    o = np.argsort(-f)[:top * 5]; P = th[o]
    f = Fg(P); lr = np.full(len(P), 0.01)
    for _ in range(150):  # backtracking ascent with central-difference gradient (3 coords)
        g = np.zeros_like(P)
        for j in idx:
            e = np.zeros(6); e[j] = 1e-6; g[:, j] = (Fg(P + e) - Fg(P - e)) / 2e-6
        nP = P + lr[:, None] * g; fn = Fg(nP); ok = fn > f
        P[ok] = nP[ok]; f[ok] = fn[ok]; lr = np.where(ok, lr * 1.3, lr * 0.4)
    o = np.argsort(-f); P = P[o]; f = f[o]
    # dedupe
    keep = []
    for i in range(len(P)):
        if all(np.abs(((P[i] - P[j] + math.pi) % (2 * math.pi)) - math.pi).max() > 0.05 for j in keep): keep.append(i)
        if len(keep) >= top: break
    return P[keep][:, idx], f[keep]
PA, fA = group_max([0, 1, 2]); PB, fB = group_max([3, 4, 5])
print(f'group maxima: A top {fA[:3]} ({len(PA)} kept, last {fA[-1]:.4f}); B top {fB[:3]} ({len(PB)} kept, last {fB[-1]:.4f}) [{time.time()-t0:.0f}s]', flush=True)
prod = np.array([np.r_[x, y] for x in PA[:60] for y in PB[:60]])
def truemax(a, extra):
    cand = np.vstack([prod, rs.uniform(-math.pi, math.pi, (NR, 6)), extra])
    f = F(cand, a); o = np.argsort(-f)[:NTOP]; P, f = refine(cand[o], a)
    o = np.argsort(-f); return float(f[o[0]]), P[o]
s_true, P = truemax(a, np.zeros((0, 6)))
print(f'init true max {s_true:.6f} (separable cert 5.189) [{time.time()-t0:.0f}s]', flush=True)
pts = np.vstack([prod, P, rs.uniform(-math.pi, math.pi, (3000, 6))]); pts = pts[np.argsort(-F(pts, a))[:PMAX]]
best_a, best_s = a.copy(), s_true; r = float(os.environ.get('R0', '0.01')); hist = []
for it in range(400):
    C = np.cos(pts @ Kf.T); Bp = Bv(pts)
    res = linprog(np.r_[np.zeros(len(K)), 1.0], A_ub=np.hstack([-C, -np.ones((len(pts), 1))]), b_ub=-Bp,
                  bounds=[(x - r, x + r) for x in best_a] + [(None, None)], method='highs')
    if res.status != 0: print('LP status', res.status, res.message, flush=True); r *= 0.5; continue
    an = res.x[:-1]; s_lp = float(res.x[-1]); s_new, Pn = truemax(an, P[:1500])
    acc = s_new < best_s
    if acc: best_a, best_s = an, s_new; r = min(r * 1.5, 1.0)
    else: r *= 0.5
    hist.append((s_lp, s_new, r))
    print(f' it {it}: LP {s_lp:.6f} true {s_new:.6f} best {best_s:.6f} {"ACC" if acc else "rej"} r={r:.4g} pts={len(pts)} [{time.time()-t0:.0f}s]', flush=True)
    Fp = Bp - C @ best_a; keepm = Fp > best_s - 0.4
    pts = np.vstack([pts[keepm], Pn[:1200], P[:300]])
    if len(pts) > PMAX: pts = pts[np.argsort(-F(pts, best_a))[:PMAX]]
    if acc: P = Pn
    cx = len(K0); out = dict(L='7/5', EPS=EPS, DC=DC, primes=ps, vis=[(int(i), int(k), float(c)) for (i, k, c) in vis], K=K.tolist(),
        a_lp=best_a.tolist(), A_lp=s_lp, A_true=best_s, n_cross=len(cross), cross_l1=float(np.abs(best_a[cx:]).sum()),
        sep_moved_l1=float(np.abs(best_a[:cx] - a0).sum()), hist=hist, label='NUMERICAL (6-D sup estimate; not B&B)')
    json.dump(out, open(f'joint_{tag}.json', 'w'))
    if r < 1e-4 or time.time() - t0 > budget: break
print('DONE best', best_s, 'T(beta=0.15) =', math.ceil(2 * math.pi * math.exp(best_s + 0.15) + 1), flush=True)
