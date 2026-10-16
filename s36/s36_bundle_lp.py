"""Grok S36 H36b: proximal (trust-region) cutting-plane LP for a 3-prime torus lift (fixes the S34/S35
exchange-LP oscillation at {2,3,5} D=8, L=1.3). True max estimated by a full 3D grid + gradient refinement.
NUMERICAL design only; rigour via s35_interval_bb.py. RH OPEN.
usage: s36_bundle_lp.py L D budget_s tag primes [init.json]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's28'))
from s28_torus_lp import setup, chars, Bval, Fval, Fgrad
from scipy.optimize import linprog
L = float(sys.argv[1]); D = int(sys.argv[2]); budget = float(sys.argv[3]); tag = sys.argv[4]
plist = [int(x) for x in sys.argv[5].split(',')]
init = sys.argv[6] if len(sys.argv) > 6 else None
rs = np.random.default_rng(1); t0 = time.time()
ps, lp, vis = setup(L, plist); n = len(ps)
K = chars(lp, L, D, 6 * (max(k for (_, k, _) in vis) + 1))
print(f'L={L} D={D} primes={ps} chars={len(K)}', flush=True)
a = np.zeros(len(K))
if init:
    d = list(json.load(open(init)).values())[0]
    Ki = [tuple(k) for k in d['K']]; m = {k: x for k, x in zip(Ki, d['a_lp'])}
    a = np.array([m.get(tuple(k), 0.0) for k in K.tolist()])
NG = int(os.environ.get('NGRID', '96'))
g1 = (np.arange(NG) + 0.5) * 2 * math.pi / NG - math.pi
def refine(th, a, iters=60):
    lr = 0.01
    f = Fval(th, vis, K, a)
    for _ in range(iters):
        g = Fgrad(th, vis, K, a); nt = th + lr * g
        fn = Fval(nt, vis, K, a); ok = fn > f
        th[ok] = nt[ok]; f[ok] = fn[ok]
    return th, f
def truemax(a, keep=3000):
    best = []
    # grid over theta_1 in [0,pi] (F even), full in others
    for i in range(NG // 2, NG):
        G = np.stack(np.meshgrid([g1[i]], g1, g1, indexing='ij'), -1).reshape(-1, 3)
        f = Fval(G, vis, K, a); idx = np.argsort(-f)[:200]; best.append((G[idx], f[idx]))
    P = np.vstack([b[0] for b in best]); f = np.concatenate([b[1] for b in best])
    idx = np.argsort(-f)[:keep]; P = P[idx].copy()
    P, f = refine(P, a)
    return float(f.max()), P[np.argsort(-f)[:keep]]
s_true, P = truemax(a)
pts = np.vstack([rs.uniform(-math.pi, math.pi, (6000, n)), P])
best_a, best_s = a.copy(), s_true; r = 0.05; hist = []
print(f'init true max {s_true:.6f} [{time.time()-t0:.0f}s]', flush=True)
for it in range(200):
    C = np.cos(pts @ K.T); Bp = Bval(pts, vis)
    res = linprog(np.r_[np.zeros(len(K)), 1.0], A_ub=np.hstack([-C, -np.ones((len(pts), 1))]), b_ub=-Bp,
                  bounds=[(x - r, x + r) for x in best_a] + [(None, None)], method='highs')
    an = res.x[:-1]; s_lp = float(res.x[-1])
    s_new, Pn = truemax(an)
    acc = s_new < best_s
    if acc: best_a, best_s = an, s_new; r = min(r * 1.5, 1.0)
    else: r *= 0.5
    hist.append((s_lp, s_new, r))
    print(f' it {it}: LP {s_lp:.6f} true {s_new:.6f} best {best_s:.6f} {"ACC" if acc else "rej"} r={r:.4g} pts={len(pts)} [{time.time()-t0:.0f}s]', flush=True)
    # keep only near-active points to keep LP small
    Fp = Bp - C @ best_a; keepm = Fp > best_s - 0.6
    pts = np.vstack([pts[keepm], Pn[:1500], P[:500] if acc else np.zeros((0, n))])
    P = Pn if acc else P
    out = dict(L=L, D=D, primes=ps, vis=[(int(i), int(k), float(c)) for (i, k, c) in vis], K=K.tolist(), a_lp=best_a.tolist(),
               A_lp=s_lp, A_true=best_s, hist=hist, NG=NG)
    json.dump({f'L{L}_D{D}': out}, open(f'bundle_{tag}.json', 'w'))
    if r < 1e-4 or time.time() - t0 > budget: break
print('DONE best', best_s, flush=True)
