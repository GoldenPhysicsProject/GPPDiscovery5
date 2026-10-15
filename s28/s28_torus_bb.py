"""Grok S28: rigorous upper bound for A = sup_{theta in T^n} F(theta),
F = B - m,  B = sum c_{p,k} cos(k theta_p) (visible comb),  m = sum_j a_j cos(K_j . theta) (invisible lift).
Branch and bound over boxes (centre c, half-widths d):
  F(x) <= F(c) + sum_i |dF/dx_i(c)| d_i + (1/2) sum_{i,l} H_il d_i d_l,   H_il = sum_terms |coef| |k_i| |k_l|
(Taylor with Lagrange remainder; |d^2 cos| <= 1).  Symmetry F(-x) = F(x) => x_1 in [0, pi].
Floating point: float64 evaluation; an explicit a-priori rounding bound
  E_round = 8 eps * sum_terms |coef| (1 + |k|_1 pi + n_terms) + 1e-12
is added to every box bound (assumes libm cos/sin within 2 ulp).  Coefficients are first rounded to
decimal strings (a_q), which are exactly what the Arb certificate uses.  RH OPEN."""
import numpy as np, json, sys, math, time

def load(fn, digits=9):
    d = json.load(open(fn)); ps = d['primes']; n = len(ps)
    rows = []; coefs = []
    for (i, k, c) in d['vis']:
        kv = [0] * n; kv[i] = k; rows.append(kv); coefs.append(c)
    a_q = ['%.*e' % (digits, x) for x in d['a']]
    for kv, x in zip(d['K'], a_q):
        rows.append(kv); coefs.append(-float(x))
    return d, ps, np.array(rows, float), np.array(coefs), a_q

def bb(Kmat, coef, n, target, d0=math.pi / 8, maxboxes=6_000_000, verbose=True):
    eps = 2.0 ** -52
    absc = np.abs(coef); k1 = np.abs(Kmat).sum(1)
    Eround = 8 * eps * float((absc * (1 + k1 * math.pi + len(coef))).sum()) + 1e-12
    Kabs = np.abs(Kmat); H = (Kabs.T * absc) @ Kabs
    m1 = int(round(math.pi / (2 * d0))); mo = 2 * m1
    g1 = (np.arange(m1) + 0.5) * (math.pi / m1); go = -math.pi + (np.arange(mo) + 0.5) * (2 * math.pi / mo)
    C = np.array(np.meshgrid(*([g1] + [go] * (n - 1)), indexing='ij')).reshape(n, -1).T
    D = np.full((len(C), n), math.pi / (2 * m1))
    best = -1e9; nproc = 0; it = 0
    while len(C):
        it += 1
        ub = np.empty(len(C)); gabs = np.empty_like(C)
        for s in range(0, len(C), 200000):
            ph = C[s:s+200000] @ Kmat.T
            F = np.cos(ph) @ coef
            grad = -(np.sin(ph) * coef) @ Kmat
            Dd = D[s:s+200000]
            best = max(best, F.max())
            ub[s:s+200000] = F + (np.abs(grad) * Dd).sum(1) + 0.5 * np.einsum('bi,il,bl->b', Dd, H, Dd) + Eround
            gabs[s:s+200000] = np.abs(grad)
        nproc += len(C)
        keep = ub > target
        if verbose: print(f'  it {it}: boxes {len(C)} keep {keep.sum()} best F {best:.6f} max ub {ub.max():.6f} min d {D.min():.2e}', flush=True)
        if best > target: return None, best, nproc, Eround
        C = C[keep]; D = D[keep]; gabs = gabs[keep]
        if len(C) == 0: break
        if len(C) * 2 > maxboxes: return None, best, nproc, Eround
        j = np.argmax(D * (gabs + D @ H), axis=1)
        D2 = D.copy(); D2[np.arange(len(D)), j] /= 2
        off = np.zeros_like(D); off[np.arange(len(D)), j] = D2[np.arange(len(D)), j]
        C = np.vstack([C - off, C + off]); D = np.vstack([D2, D2])
    return target, best, nproc, Eround

if __name__ == '__main__':
    fn = sys.argv[1]; margin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.005
    d, ps, Kmat, coef, a_q = load(fn); n = len(ps)
    target = float('%.4f' % (d['A_max_est'] + margin)); t0 = time.time()
    while True:
        r, best, nproc, Er = bb(Kmat, coef, n, target)
        if r is not None: break
        target = float('%.4f' % (max(target, best) + margin)); print('raise target', target, flush=True)
    print(f'{fn}: sup_torus F <= {target} (best sample {best:.6f}); boxes processed {nproc}; Eround {Er:.2e} [{time.time()-t0:.0f}s]')
    out = dict(L=d['L'], primes=ps, vis=d['vis'], K=d['K'], a_q=a_q, A_cert=repr(target), A_best_sample=best,
               method='float64 B&B, 2nd-order Taylor, explicit rounding margin', boxes=nproc, Eround=Er)
    json.dump(out, open(fn.replace('.json', '_cert.json'), 'w'))
