"""Grok S29 H29a (corrected): eigenfunction-seeded torus LP for invisible lift.
DEAD END logged: plain Chebyshev LP on a finite t-interval overfits (A_grid~3.07 at L=1.19
but true almost-periodic max ~14). Proper route = Bohr/torus cutting planes (S28), seeded by
the dual measure |F(t)|^2 of the top eigenfunction of P_L: theta_p = (t log p) mod 2pi at
weight peaks. NUMERICAL. RH OPEN."""
import numpy as np, math, json, sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's28'))
from s28_pmax import Pmat, mangoldt_terms, phi
from s28_torus_lp import setup, chars, Bval, Fval, Fgrad, local_max
from scipy.optimize import linprog
from numpy.polynomial.legendre import leggauss

def eigen_torus_seeds(L, N=100, n_seeds=8000, t_max=200.0, nt=20000):
    P = Pmat(L, N)
    ev, V = np.linalg.eigh(P)
    lam = float(ev[-1]); v = V[:, -1]
    ps, lp, vis = setup(L)
    n = len(ps)
    # synthesize F
    xg, wg = leggauss(max(4*N, 200))
    x = L*xg; w = L*wg
    f = sum(v[i]*phi(i, x, L) for i in range(len(v)))
    t = np.linspace(0, t_max, nt)
    F = np.exp(-1j*np.outer(t, x)) @ (f*w)
    wgt = np.abs(F)**2
    # top weight indices
    idx = np.argsort(wgt)[-n_seeds:]
    t_star = t[idx]
    # map to torus
    th = np.mod(np.outer(t_star, lp) + math.pi, 2*math.pi) - math.pi  # in [-pi,pi]
    # also add random for coverage
    return lam, ps, lp, vis, th, float(ev[-1])

def solve_seeded(L, D=3, N=100, npts0=15000, rounds=20, amax=20.0, seed=1, gap=3e-4):
    rs = np.random.default_rng(seed)
    lam, ps, lp, vis, th_seed, _ = eigen_torus_seeds(L, N=N, n_seeds=6000)
    n = len(ps)
    kmax_axis = 6 * (max(k for (_, k, _) in vis) + 1)
    K = chars(lp, L, D, kmax_axis)
    pts = np.vstack([th_seed, rs.uniform(-math.pi, math.pi, (npts0, n))])
    a = np.zeros(len(K)); hist = []; smax = 1e9; s = 1e9
    print(f'L={L} D={D} primes={ps} chars={len(K)} seeds={len(th_seed)} lam_max={lam:.6f}', flush=True)
    for r in range(rounds):
        C = np.cos(pts @ K.T); Bp = Bval(pts, vis)
        Aub = np.hstack([-C, -np.ones((len(pts), 1))])
        res = linprog(np.r_[np.zeros(len(K)), 1.0], A_ub=Aub, b_ub=-Bp,
                      bounds=[(-amax, amax)] * len(K) + [(None, None)], method='highs')
        if not res.success:
            print('LP fail', res.message); break
        a = res.x[:-1]; s = float(res.x[-1])
        cand = rs.uniform(-math.pi, math.pi, (50000, n))
        # also re-seed from eigen peaks each round (perturb)
        cand = np.vstack([cand, th_seed + rs.normal(0, 0.1, th_seed.shape),
                          pts[np.argsort(-(Bp - C@a))[:3000]]])
        cand = local_max(cand, vis, K, a)
        Fc = Fval(cand, vis, K, a); smax = float(Fc.max())
        hist.append((s, smax))
        print(f'  round {r}: LP s={s:.6f} true max>={smax:.6f} pts={len(pts)}', flush=True)
        if smax - s < gap: break
        pts = np.vstack([pts, cand[np.argsort(-Fc)[:4000]]])
    return dict(L=L, D=D, primes=ps, n_chars=len(K), lam_max=lam,
                A_lp=s, A_max_est=smax, hist=hist,
                K=K.tolist(), a=a.tolist(),
                vis=[(int(i),int(k),float(c)) for (i,k,c) in vis],
                per_prime_ref={1.0:3.047608, 1.19:3.919408}.get(round(L,2)),
                torus_S28_ref={1.0:2.95, 1.19:3.602}.get(round(L,2)),
                T1=2*math.pi*math.exp(smax),
                T1_floor=2*math.pi*math.exp(lam),
                T1_per_prime=2*math.pi*math.exp({1.0:3.047608, 1.19:3.919408}.get(round(L,2), smax)))

if __name__ == '__main__':
    os.environ['OMP_NUM_THREADS'] = '1'
    out = {}
    t0 = time.time()
    for L, D in [(1.0, 3), (1.19, 3)]:
        key = f'L{L}_D{D}'
        out[key] = solve_seeded(L, D=D, N=80 if L<=1 else 100, rounds=18, gap=4e-4)
        print(f"RESULT {key}: A~{out[key]['A_max_est']:.5f} (LP {out[key]['A_lp']:.5f}) "
              f"lam={out[key]['lam_max']:.5f} S28torus={out[key]['torus_S28_ref']} "
              f"per_prime={out[key]['per_prime_ref']} T1={out[key]['T1']:.1f}", flush=True)
    json.dump(out, open('h29a_krein_torus.json','w'))
    # also write a torus-cert-shaped JSON for L=1.19 if improved
    best = out.get('L1.19_D3')
    if best and best['A_max_est'] < 3.55:
        cert = dict(L=1.19, D=3, A_lp=best['A_lp'], A_max_est=best['A_max_est'],
                    A_cert=best['A_max_est']+5e-4,  # margin pending B&B
                    primes=best['primes'], K=best['K'], a=best['a'], vis=best['vis'],
                    note='NUMERICAL torus max; not B&B certified. Eigenfunction-seeded S29.')
        json.dump(cert, open('torus_L1.19_D3_eigen.json','w'))
        print('Wrote torus_L1.19_D3_eigen.json (NUMERICAL, not B&B)')
    print(f'done in {time.time()-t0:.0f}s')
