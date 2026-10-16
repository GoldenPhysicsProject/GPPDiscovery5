"""Tighten L=1 D=4 Fejér–Riesz torus A with more cutting-plane rounds + larger Bohr audit.
NUMERICAL. RH OPEN."""
import numpy as np, math, json, sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's28'))
from s28_torus_lp import setup, chars, Bval, Fval, local_max
from scipy.optimize import linprog

def refine(path='h30a_fejer_riesz.json', key='L1.0_D4', rounds=20, gap=2e-4, n_rand=400000):
    data = json.load(open(path))
    r = data[key]
    L, D = r['L'], r['D']
    ps, lp, vis = setup(L)
    K = np.array(r['K'])
    a = np.array(r['a_lp'])
    rs = np.random.default_rng(7)
    n = K.shape[1]
    # seed from previous + random
    pts = rs.uniform(-math.pi, math.pi, (30000, n))
    # add peaks of current F
    cand0 = rs.uniform(-math.pi, math.pi, (80000, n))
    Fc0 = Fval(cand0, vis, K, a)
    pts = np.vstack([pts, cand0[np.argsort(-Fc0)[:8000]]])
    hist = []
    s = smax = 1e9
    for rnd in range(rounds):
        C = np.cos(pts @ K.T); Bp = Bval(pts, vis)
        Aub = np.hstack([-C, -np.ones((len(pts), 1))])
        res = linprog(np.r_[np.zeros(len(K)), 1.0], A_ub=Aub, b_ub=-Bp,
                      bounds=[(-25, 25)] * len(K) + [(None, None)], method='highs')
        if not res.success:
            print('LP fail', res.message); break
        a = res.x[:-1]; s = float(res.x[-1])
        cand = rs.uniform(-math.pi, math.pi, (100000, n))
        cand = np.vstack([cand, pts[np.argsort(-(Bp - C @ a))[:4000]]])
        cand = local_max(cand, vis, K, a, iters=300, lr=0.025)
        Fc = Fval(cand, vis, K, a); smax = float(Fc.max())
        hist.append((s, smax))
        print(f'  refine {rnd}: LP={s:.6f} true>={smax:.6f} pts={len(pts)}', flush=True)
        if smax - s < gap: break
        pts = np.vstack([pts, cand[np.argsort(-Fc)[:5000]]])
    # independent Bohr audit
    audit = rs.uniform(-math.pi, math.pi, (n_rand, n))
    Fa = Fval(audit, vis, K, a)
    top = audit[np.argsort(-Fa)[:15000]]
    top = local_max(top, vis, K, a, iters=400)
    A_true = float(max(Fa.max(), Fval(top, vis, K, a).max()))
    out = dict(L=L, D=D, A_lp=s, A_max_est=smax, A_true=A_true, hist=hist,
               lam_max=r['lam_max'], per_prime=r['per_prime'], s28_torus=r['s28_torus'],
               ratio_to_floor=A_true / r['lam_max'],
               gap_to_floor=A_true - r['lam_max'],
               K=K.tolist(), a=a.tolist(),
               vis=[(int(i),int(k),float(c)) for (i,k,c) in vis],
               primes=ps,
               note='NUMERICAL refined D=4; not B&B certified.')
    json.dump(out, open('h30a_L1_D4_refined.json', 'w'), indent=1)
    print(f'REFINED A_true={A_true:.6f} LP={s:.6f} ratio={out["ratio_to_floor"]:.4f}')
    return out

if __name__ == '__main__':
    os.environ['OMP_NUM_THREADS']='1'
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    refine()
