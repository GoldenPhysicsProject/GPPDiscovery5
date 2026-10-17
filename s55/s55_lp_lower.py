"""Grok S55 H55a: NUMERICAL lower bound for the best separable torus-lift constant over a fixed character set.
For a finite point set P on the torus, min_a max_{theta in P} [B(theta) - sum_k a_k cos(k.theta)] <= min_a sup_theta (...)
= A*(char set). So the LP value is a lower bound on the best achievable A for that part (up to LP tolerance).
Characters: |k|_1 <= D cross terms and axis terms, all with |k.logp| >= 2L + EPS (buffer rule eps*T >~ 5).
usage: s55_lp_lower.py L D EPS primes NPTS [seed]   RH OPEN."""
import sys, os, math, json, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's28'))
from s28_torus_lp import setup, chars, Bval
from scipy.optimize import linprog
L = float(sys.argv[1]); D = int(sys.argv[2]); EPS = float(sys.argv[3]); plist = [int(x) for x in sys.argv[4].split(',')]
NP = int(sys.argv[5]); seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0
t0 = time.time(); rs = np.random.default_rng(seed)
ps, lp, vis = setup(L, plist); n = len(ps)
K = chars(lp, L + EPS / 2, D, 6 * (max(k for (_, k, _) in vis) + 1))
pts = rs.uniform(-math.pi, math.pi, (NP, n))
# add the points where B is large (B is maximized at theta=0 etc.): include the lattice of small rational points
g = np.linspace(-math.pi, math.pi, 25)[:-1]
G = np.stack(np.meshgrid(*([g] * n), indexing='ij'), -1).reshape(-1, n)
pts = np.vstack([pts, G])
for it in range(3):
    C = np.cos(pts @ K.T); Bp = Bval(pts, vis)
    res = linprog(np.r_[np.zeros(len(K)), 1.0], A_ub=np.hstack([-C, -np.ones((len(pts), 1))]), b_ub=-Bp,
                  bounds=[(None, None)] * len(K) + [(None, None)], method='highs')
    a = res.x[:-1]; lb = float(res.x[-1])
    # evaluate achieved sup on a fresh dense random sample (upper estimate of A for this a, still only a sample)
    Q = rs.uniform(-math.pi, math.pi, (200000, n)); F = Bval(Q, vis) - np.cos(Q @ K.T) @ a
    j = np.argsort(-F)[:3000]
    print(json.dumps(dict(L=L, D=D, EPS=EPS, primes=ps, nchars=len(K), minfreq=float(np.min(np.abs(K @ lp))), npts=len(pts),
                          LP_lower=lb, sample_sup=float(F.max()), it=it, t=round(time.time() - t0))), flush=True)
    pts = np.vstack([pts, Q[j]])
