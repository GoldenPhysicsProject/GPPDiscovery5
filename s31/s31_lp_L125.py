"""Grok S31 H31b: D=4 cross-prime torus LP at L=1.25 (primes 2,3,5,7,11), time-capped, saves each round.
Same model as s28_torus_lp.solve. NUMERICAL design; rigour via s31_bb_D4-style B&B.  RH OPEN."""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's28'))
from s28_torus_lp import setup, chars, Bval, Fval, local_max
from scipy.optimize import linprog
L = float(sys.argv[1]); D = int(sys.argv[2]); budget = float(sys.argv[3]); tag = sys.argv[4]
rs = np.random.default_rng(0); t0 = time.time()
plist = [int(x) for x in sys.argv[5].split(',')] if len(sys.argv) > 5 else None
ps, lp, vis = setup(L, plist); n = len(ps)
K = chars(lp, L, D, 6 * (max(k for (_, k, _) in vis) + 1))
print(f'L={L} D={D} primes={ps} vis={len(vis)} chars={len(K)} per_prime_raw={sum(c for *_, c in vis):.4f}', flush=True)
pts = rs.uniform(-math.pi, math.pi, (15000, n)); hist = []
for r in range(40):
    C = np.cos(pts @ K.T); Bp = Bval(pts, vis)
    res = linprog(np.r_[np.zeros(len(K)), 1.0], A_ub=np.hstack([-C, -np.ones((len(pts), 1))]), b_ub=-Bp,
                  bounds=[(-20, 20)] * len(K) + [(None, None)], method='highs')
    a = res.x[:-1]; s = float(res.x[-1])
    cand = np.vstack([rs.uniform(-math.pi, math.pi, (30000, n)), pts[np.argsort(-Bp + C @ a)[:3000]] + rs.normal(0, 0.05, (3000, n))])
    cand = local_max(cand, vis, K, a, iters=150)
    Fc = Fval(cand, vis, K, a); smax = float(Fc.max()); hist.append((s, smax))
    print(f'  round {r}: LP s={s:.6f} true>={smax:.6f} pts={len(pts)} [{time.time()-t0:.0f}s]', flush=True)
    out = dict(L=L, D=D, primes=ps, vis=[(int(i), int(k), float(c)) for (i, k, c) in vis], K=K.tolist(), a_lp=a.tolist(),
               A_lp=s, A_true=smax, hist=hist, per_prime=float(sum(c for *_, c in vis)), lam_max=None)
    json.dump({f'L{L}_D{D}': out}, open(f'lp_{tag}.json', 'w'))
    if smax - s < 2e-3 or time.time() - t0 > budget: break
    pts = np.vstack([pts, cand[np.argsort(-Fc)[:3000]]])
print('DONE', flush=True)
