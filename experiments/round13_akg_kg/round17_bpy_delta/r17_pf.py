"""Round 17, part A2 (consecutive points kept >= 0.02 apart to avoid cancellation artifacts): targeted PF_k search (k=2..5) for the BPY log-density of ZETA, EPS, DH.
Minimises the normalised minor det[Phi(x_i - y_j)] / prod Phi(x_i - y_i) over increasing x, y (random starts +
coordinate descent). Continuous evaluation of Phi (theta series). NUMERICAL. Output r17_pf.json"""
import sys, os, json, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from laws17 import Model
mp.mp.dps = 45
out = {}
def run(name, k, starts=60, iters=150):
    M = Model(name); cache = {}
    def P(t):
        t = abs(mp.mpf(t))
        if t not in cache: cache[t] = M.phi(mp.mpf(t))
        return cache[t]
    def f(x, y):
        if any(x[i+1] - x[i] < 0.02 for i in range(k-1)) or any(y[i+1] - y[i] < 0.02 for i in range(k-1)): return mp.inf
        d = mp.det(mp.matrix([[P(a - b) for b in y] for a in x]))
        sc = mp.fprod([P(x[i] - y[i]) for i in range(k)])
        return d/abs(sc)
    rng = random.Random(170 + k); best = None
    for _ in range(starts):
        w = rng.choice([0.3, 0.6, 1.0, 1.5])
        x = [rng.uniform(-w, w) + 0.05*i for i in range(k)]; x.sort(); x = [x[i] + 0.03*i for i in range(k)]
        y = [rng.uniform(-w, w) for i in range(k)]; y.sort(); y = [y[i] + 0.03*i for i in range(k)]
        v = f(x, y); step = w/4
        for _ in range(iters):
            imp = False
            for arr in (x, y):
                for i in range(k):
                    for sgn in (1, -1):
                        old = arr[i]; arr[i] = old + sgn*step; nv = f(x, y)
                        if nv < v: v = nv; imp = True
                        else: arr[i] = old
            if not imp:
                step /= 2
                if step < 1e-4: break
        if best is None or v < best[0]: best = (v, list(x), list(y))
    return dict(min_normalised_minor=mp.nstr(best[0], 6), x=[round(a, 5) for a in best[1]], y=[round(b, 5) for b in best[2]],
                verdict=('FAIL' if best[0] < -mp.mpf(10)**-25 else 'pass (no negative minor found)'))
for name in ['ZETA', 'DH', 'EPS']:
    out[name] = {}
    for k in [2, 3, 4, 5]:
        r = run(name, k, starts=(40 if k < 5 else 30)); out[name][f'PF{k}'] = r
        print(name, k, r, flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r17_pf.json'), 'w'), indent=1)
