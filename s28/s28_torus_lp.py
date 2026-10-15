"""Grok S28 H28a: cross-prime invisible-frequency lift (torus LP).
supp f in [-L,L]  =>  cos(u t) is invisible in Q for |u| >= 2L.  Write theta_p = t log p.
Visible comb B(theta) = sum_p sum_{k log p < 2L} c_{p,k} cos(k theta_p),  c_{p,k} = 2 log p / p^{k/2}.
Add m(theta) = sum_{k in I} a_k cos(k.theta), k in Z^n with |k.log p| >= 2L (invisible).
Then comb(t) - m(t) <= A := sup_{torus} (B - m)  (sup over t <= sup over torus).
Primal LP (cutting planes on the torus): minimise A over a_k.  Dual: probability measure on the
torus with Fourier support in the visible set (B only sees visible frequencies).
Per-prime lift (S27) = special case with k on the axes.  Float design only; rigour in s28_torus_bb.py.
RH OPEN."""
import numpy as np, math, json, sys, itertools, time, os
GAP = float(os.environ.get("TORUS_GAP", "1e-5"))
from scipy.optimize import linprog

def primes_upto(x):
    return [q for q in range(2, int(x) + 1) if all(q % r for r in range(2, int(q**.5) + 1))]

def setup(L, plist=None):
    ps = [p for p in primes_upto(math.exp(2 * L) + 1) if math.log(p) < 2 * L]
    if plist is not None: ps = plist
    lp = np.array([math.log(p) for p in ps])
    vis = []  # (prime index, k, c)
    for i, p in enumerate(ps):
        k = 1
        while k * lp[i] < 2 * L:
            vis.append((i, k, 2 * lp[i] / p ** (k / 2))); k += 1
    return ps, lp, vis

def chars(lp, L, D, kmax_axis):
    """invisible characters k (one per +-pair, first nonzero entry > 0) with |k|_1 <= D (cross) or axis k <= kmax."""
    n = len(lp); out = []
    rng = range(-D, D + 1)
    for i in range(n):
        for m in range(1, kmax_axis + 1):
            k = np.zeros(n, dtype=int); k[i] = m
            if m * lp[i] >= 2 * L * (1 + 1e-12): out.append(k)
    for k in itertools.product(rng, repeat=n):
        k = np.array(k); s1 = np.abs(k).sum()
        if (k != 0).sum() < 2 or s1 > D: continue
        nz = k[k != 0]
        if nz[0] < 0: continue
        if abs(k @ lp) >= 2 * L * (1 + 1e-12): out.append(k)
    return np.array(out)

def Bval(th, vis):
    return sum(c * np.cos(k * th[:, i]) for (i, k, c) in vis)

def Fval(th, vis, K, a):
    return Bval(th, vis) - np.cos(th @ K.T) @ a

def Fgrad(th, vis, K, a):
    g = np.zeros_like(th)
    for (i, k, c) in vis: g[:, i] -= c * k * np.sin(k * th[:, i])
    g += (np.sin(th @ K.T) * a) @ K
    return g

def local_max(th, vis, K, a, iters=200, lr=0.02):
    # projected-free gradient ascent with backtracking (Adam-like simple)
    for it in range(iters):
        g = Fgrad(th, vis, K, a)
        th = th + lr * g / (1 + np.abs(g).max(axis=1, keepdims=True) * 0.0 + 1)
    return th

def solve(L, D=3, kmax_axis=None, plist=None, npts=20000, rounds=25, amax=20.0, seed=0, verbose=True):
    rs = np.random.default_rng(seed)
    ps, lp, vis = setup(L, plist)
    n = len(ps)
    if kmax_axis is None: kmax_axis = 6 * (max(k for (_, k, _) in vis) + 1)
    K = chars(lp, L, D, kmax_axis)
    pts = rs.uniform(-math.pi, math.pi, (npts, n))
    a = np.zeros(len(K)); hist = []
    for r in range(rounds):
        C = np.cos(pts @ K.T); Bp = Bval(pts, vis)
        Aub = np.hstack([-C, -np.ones((len(pts), 1))])
        res = linprog(np.r_[np.zeros(len(K)), 1.0], A_ub=Aub, b_ub=-Bp,
                      bounds=[(-amax, amax)] * len(K) + [(None, None)], method='highs')
        a = res.x[:-1]; s = res.x[-1]
        # search for violations
        cand = rs.uniform(-math.pi, math.pi, (40000, n))
        cand = np.vstack([cand, pts[np.argsort(-Bp + C @ a)[:2000]] + rs.normal(0, 0.05, (2000, n))])
        cand = local_max(cand, vis, K, a)
        Fc = Fval(cand, vis, K, a); smax = Fc.max()
        hist.append((s, smax))
        if verbose: print(f'  round {r}: LP s={s:.6f} true max>={smax:.6f} chars={len(K)} pts={len(pts)}', flush=True)
        if smax - s < GAP: break
        top = cand[np.argsort(-Fc)[:3000]]
        pts = np.vstack([pts, top])
    return dict(L=L, primes=ps, vis=[(int(i), int(k), float(c)) for (i, k, c) in vis], K=K.tolist(), a=a.tolist(),
                A_lp=float(s), A_max_est=float(smax), per_prime_raw=float(sum(c for (_, _, c) in vis)))

if __name__ == '__main__':
    L = float(eval(sys.argv[1])); D = int(sys.argv[2]); tag = sys.argv[3] if len(sys.argv) > 3 else ''
    plist = None
    if len(sys.argv) > 4: plist = [int(x) for x in sys.argv[4].split(',')]
    t0 = time.time()
    out = solve(L, D=D, plist=plist)
    print(f'L={L} D={D} primes={out["primes"]} A_raw={out["per_prime_raw"]:.5f} torus A~{out["A_max_est"]:.5f} (LP {out["A_lp"]:.5f}) [{time.time()-t0:.0f}s]')
    json.dump(out, open(f'torus_L{L}_D{D}{tag}.json', 'w'))
