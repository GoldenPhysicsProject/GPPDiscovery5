"""Grok S53: = s36_onestroke_schurboost with --parities option (assemble only the requested parity; identical math).
Grok S33: = S32 schurboost + Bessel precision-boost fix (SPHJ_BOOST, default 1.75 instead of 1.5).
Grok S32: S28 one-stroke + Schur-precision boost (assemble at --prec, PD-check at --schur_prec).
Grok S28 original: fast rigorous lifted one-stroke certificate (Zhu arXiv:2608.24827 Thm 1.1 + S27/S28 invisible lift).
Same mathematics and error budget as GPPDiscovery5/s27/s27_onestroke_arb.py, re-engineered for speed:
 * assembly parallelised over panels (multiprocessing, fork); Bessel j_0..j_D computed ONCE per node for both parities;
 * partial Gram sums returned as Arb balls serialised with radius (printed ball contains the original ball);
 * PD check = blocked Schur-complement elimination in ball arithmetic (arb_mat.solve + arb_mat mul, C level):
   A11 PD (ball Cholesky) and A22 - A21 A11^{-1} A12 PD  <=>  A PD; every exact matrix in the ball passes.
 * shift chosen from an inverse-power estimate of lam_min (non-rigorous guess), then ONE rigorous check.
Lift input: --lift (S27 per-prime JSON, kappa bounded in Arb on a theta grid) or --torus (S28 torus cert JSON:
 list of integer character vectors k, coefficients a_k, rigorous upper bound A_cert of sup_torus(B - m)).
Status: computer-assisted (Arb), NOT Lean.  RH OPEN."""
import sys, math, json, time, argparse, os
import multiprocessing as mp
from flint import arb, acb, arb_mat, ctx

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's27'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's36'))
from s27_onestroke_arb import comb_terms, sph_j_all, sph_i, lgdf
from s27_onestroke_arb import sph_j

# Grok S33 fix: s27 sph_j_all boosts precision by 1.5*x+64 bits; Arb bessel_j(D+1/2, x) loses ALL precision for
# x ~ 896, D ~ 1500 at that boost (ball [+/- 1e-75], backward recurrence -> j_0 = [+/- 1e394], Gram rE = inf).
# Same algorithm with boost BOOST*x + 64 (default 1.75); every value is still an Arb ball (rigour unchanged).
BOOST = float(os.environ.get('SPHJ_BOOST', '1.75'))
def sph_j_all(x, D):
    rx = x.rad(); x = arb(x.mid())
    old = ctx.prec
    ctx.prec = old + int(BOOST * float(x.mid())) + 64
    try:
        v = [None] * (D + 2)
        v[D + 1] = sph_j(D + 1, x); v[D] = sph_j(D, x)
        for d in range(D, 0, -1):
            v[d - 1] = (2 * d + 1) / x * v[d] - v[d + 1]
    finally:
        ctx.prec = old
    return [+t + arb(0, rx) for t in v]

G = {}

def ser(x): return x.str(radius=True, more=True, n=max(30, int(ctx.prec * 0.302) + 5))
def des(s): return arb(s)

def symbol_real(t):
    S = acb(arb(1) / 4, t / 2).digamma().real - G['logpi'] - G['beta']
    for (_, c, lg) in G['terms']: S -= c * (t * lg).cos()
    for (x, fr) in G['lift_terms']: S += x * (t * fr).cos()
    return S

def symbol_cplx(z):
    v = ((acb(arb(1) / 4) + acb(0, 1) * z / 2).digamma() + (acb(arb(1) / 4) - acb(0, 1) * z / 2).digamma()) / 2 - G['logpi'] - G['beta']
    for (_, c, lg) in G['terms']: v -= c * (z * lg).cos()
    for (x, fr) in G['lift_terms']: v += x * (z * fr).cos()
    return abs(v).upper()

def work_bounds(job):
    """boundary sweep (max modulus) and real-axis sup, over ix in [i0, i1)"""
    ctx.prec = G['prec']; kind, i0, i1 = job; hb = G['hb']; b = G['b']; m = arb(0)
    for ix in range(i0, i1):
        if kind == 'top':
            xc = arb(G['x0']) + (2 * ix + 1) * arb(hb)
            for yy in (b, -b):
                v = symbol_cplx(acb(arb(xc, hb * 1.001), yy)); m = v if v > m else m
        elif kind == 'side':
            yc = -b + (2 * ix + 1) * arb(hb)
            for xx in (G['x0'], G['x1']):
                v = symbol_cplx(acb(arb(xx), arb(yc, hb * 1.001))); m = v if v > m else m
        else:
            xc = (2 * ix + 1) * arb(hb)
            v = abs(symbol_real(arb(xc.mid(), hb * 1.001))).upper(); m = v if v > m else m
    return ser(arb(m))

def work_panels(job):
    ctx.prec = G['prec']; panels = job
    L, w, gl, K = G['L'], G['w'], G['gl'], G['K']
    D = 2 * K  # max degree needed: odd sector 1 + 2(K-1) = 2K-1
    Cs = [arb_mat(K, K), arb_mat(K, K)]
    cols = [[], []]; wts = []
    def flush():
        n = len(wts)
        if n == 0: return
        for par in G['pars']:
            F = arb_mat(K, n); FW = arb_mat(n, K)
            for j in range(n):
                cj = cols[par][j]; wj = wts[j]
                for i in range(K):
                    F[i, j] = cj[i]; FW[j, i] = cj[i] * wj
            Cs[par] = Cs[par] + F * FW
        cols[0].clear(); cols[1].clear(); wts.clear()
    for ip in panels:
        lo = ip * w
        for (xk, wk) in gl:
            t = lo + w / 2 * (1 + xk)
            S = symbol_real(t)
            jall = sph_j_all(L * t, D)
            for par in G['pars']:
                cols[par].append([((-1) ** (d // 2)) * G['nrm'][d] * jall[d] for d in range(par, par + 2 * K, 2)])
            wts.append(w / 2 * wk * S / arb.pi())
            if len(wts) >= G['chunk']: flush()
    flush()
    return [([[ser(Cs[par][i, j]) for j in range(K)] for i in range(K)] if par in G['pars'] else None) for par in (0, 1)]

def pd_schur(Mrows, bs=48):
    """Rigorous: True only if every symmetric matrix in the ball matrix Mrows is positive definite."""
    M = arb_mat(Mrows)
    n = M.nrows()
    try:
        while True:
            b = min(bs, n)
            A11 = arb_mat([[M[i, j] for j in range(b)] for i in range(b)])
            # ball Cholesky of A11
            R = [[A11[i, j] for j in range(b)] for i in range(b)]
            for k in range(b):
                for j in range(k): R[k][k] -= R[k][j] * R[k][j]
                if not (R[k][k] > 0): return False
                R[k][k] = R[k][k].sqrt()
                for i in range(k + 1, b):
                    acc = R[i][k]
                    for j in range(k): acc -= R[i][j] * R[k][j]
                    R[i][k] = acc / R[k][k]
            if n == b: return True
            A12 = arb_mat([[M[i, j] for j in range(b, n)] for i in range(b)])
            A21 = arb_mat([[M[i, j] for j in range(b)] for i in range(b, n)])
            A22 = arb_mat([[M[i, j] for j in range(b, n)] for i in range(b, n)])
            X = A11.solve(A12)
            M = A22 - A21 * X
            n = n - b
    except (ZeroDivisionError, ValueError):
        return False

def lam_estimate(Atil, K, iters=60):
    """non-rigorous guess of lam_min via inverse power iteration on midpoints"""
    A = arb_mat(Atil)
    try:
        Ai = A.inv()
    except ZeroDivisionError:
        return None
    Ai = Ai.mid()
    v = arb_mat([[arb(1) / (i + 1)] for i in range(K)])
    lam = None
    for _ in range(iters):
        u = Ai * v
        nu = sum((u[i, 0] ** 2 for i in range(K)), arb(0)).sqrt()
        rq = sum((u[i, 0] * v[i, 0] for i in range(K)), arb(0)) / sum((v[i, 0] ** 2 for i in range(K)), arb(0))
        lam = 1 / rq
        v = arb_mat([[(u[i, 0] / nu).mid()] for i in range(K)])
    return float(lam.mid())

def sexp(x):
    # upper bound for exp(x), safe against underflow (exp(x) < 1e-300 when x < -690.8)
    return math.exp(x) if x > -690 else 1e-300

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', default='1'); ap.add_argument('--T', default='300'); ap.add_argument('--K', type=int, default=262)
    ap.add_argument('--w', default='1/2'); ap.add_argument('--ng', type=int, default=40); ap.add_argument('--b', default='2/5')
    ap.add_argument('--prec', type=int, default=384); ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--lift', default=None); ap.add_argument('--torus', default=None); ap.add_argument('--kgrid', type=int, default=40000)
    ap.add_argument('--chunk', type=int, default=1024); ap.add_argument('--hb', type=float, default=0.002)
    ap.add_argument('--bs', type=int, default=48); ap.add_argument('--out', default=None)
    ap.add_argument('--estimate_only', action='store_true')
    ap.add_argument('--schur_prec', type=int, default=0, help='if >0, run ball Schur at this prec (else use --prec)')
    ap.add_argument('--dump_atil', default=None, help='path prefix to dump even/odd Atil mid strings')
    ap.add_argument('--parities', default='0,1', help='0=even,1=odd')
    ap.add_argument('--schur_bs', default='48,24,32,64,16', help='comma list of block sizes to try')
    a = ap.parse_args()
    ctx.prec = a.prec; T0 = time.time()
    L = arb(a.L); T = arb(a.T); w = arb(a.w); b = arb(a.b)
    Lf, Tf, wf, bf = float(L.mid()), float(T.mid()), float(w.mid()), float(b.mid())
    terms = comb_terms(Lf)
    AL = sum((c for (_, c, _) in terms), arb(0)); lift_terms = []; info = {}
    if a.lift:  # S27 per-prime lift, kappa in Arb on theta grid (as s27)
        lj = json.load(open(a.lift)); AL = arb(0); kap_d = {}
        for ps, d in lj['primes'].items():
            p = int(ps); lp = arb(p).log(); ks = d['ks']; cs = [2 * lp / arb(p) ** (arb(k) / 2) for k in ks]
            inv = d['inv']; av = [arb('%.12e' % x) for x in d['a']]
            for k in inv: assert arb(k) * lp >= 2 * L
            ng = a.kgrid; hw = arb.pi() / (2 * ng); kap = None
            for i in range(ng):
                th = arb(((2 * i + 1) * hw).mid(), (hw * 1.0001).mid())
                v = sum((c * (k * th).cos() for k, c in zip(ks, cs)), arb(0)) - sum((x * (k * th).cos() for k, x in zip(inv, av)), arb(0))
                u = v.upper(); kap = u if kap is None or u > kap else kap
            AL += arb(kap); kap_d[p] = float(kap)
            lift_terms += [(x, k * lp) for k, x in zip(inv, av)]
        info['kappas'] = kap_d
    if a.torus:  # S28 torus lift; A_cert proven by s28_torus_bb.py
        tj = json.load(open(a.torus)); ps = tj['primes']
        lps = [arb(p).log() for p in ps]
        assert [n for n, _, _ in terms if all(n % q for q in range(2, n)) ] == [p for p in ps], 'prime set mismatch'
        for k, x in zip(tj['K'], tj['a_q']):
            u = sum((kk * lp for kk, lp in zip(k, lps)), arb(0))
            assert abs(u) >= 2 * L, ('visible character!', k)
            lift_terms.append((arb(x), abs(u)))
        AL = arb(tj['A_cert'])
        info['torus'] = a.torus; info['A_cert'] = tj['A_cert']
    beta = (T / (2 * arb.pi())).log() - 1 / T - AL
    print(f'L={a.L} T={a.T} K={a.K} comb={[n for n,_,_ in terms]} A={float(AL.mid()):.6f} beta={float(beta.mid()):.6f} lift_terms={len(lift_terms)}', flush=True)
    assert T > arb(15) / 4 and beta > 0, 'beta <= 0'
    npan = int(round(Tf / wf)); assert abs(npan * wf - Tf) < 1e-12
    gl = [arb.legendre_p_root(a.ng, k, weight=True) for k in range(a.ng)]
    s = 4 * bf / wf; rho = (s + math.sqrt(s * s + 4)) / 2 * (1 - 1e-9)
    ovh = (wf / 2) * ((rho + 1 / rho) / 2 - 1)
    K = a.K
    nrm = [(2 * L * (2 * d + 1)).sqrt() for d in range(2 * K + 2)]
    G.update(dict(prec=a.prec, L=L, w=w, gl=gl, K=K, terms=terms, lift_terms=lift_terms, beta=beta, logpi=arb.pi().log(),
                  nrm=nrm, pars=tuple(int(x) for x in a.parities.split(',')), chunk=a.chunk, hb=a.hb, b=b, x0=-ovh, x1=Tf + ovh))
    pool = mp.get_context('fork').Pool(a.workers)
    # --- bounds
    t0 = time.time(); hb = a.hb
    ntop = int(math.ceil((Tf + 2 * ovh) / (2 * hb))); nside = int(math.ceil(2 * bf / (2 * hb))); nreal = int(math.ceil(Tf / (2 * hb)))
    def split(kind, n, parts):
        st = max(1, n // parts); return [(kind, i, min(n, i + st)) for i in range(0, n, st)]
    rs = pool.map(work_bounds, split('top', ntop, 4 * a.workers) + split('side', nside, 1))
    MS = max((des(x).upper() for x in rs), key=lambda z: float(z)); MS = arb(MS)
    rr = pool.map(work_bounds, split('real', nreal, 4 * a.workers))
    M0 = arb(max((des(x).upper() for x in rr), key=lambda z: float(z)))
    print(f'M_S <= {float(MS):.4f}  M0 <= {float(M0):.4f} [{time.time()-t0:.0f}s]', flush=True)
    # --- assembly
    t0 = time.time()
    nj = a.workers
    jobs = [list(range(i, npan, nj)) for i in range(nj)]  # interleaved panels for load balance
    Csum = [[[arb(0)] * K for _ in range(K)] for _ in (0, 1)]
    for part in pool.imap_unordered(work_panels, jobs):
        for par in G['pars']:
            for i in range(K):
                row = Csum[par][i]; pr = part[par][i]
                for j in range(K): row[j] += des(pr[j])
    pool.close()
    print(f'assembly done: {npan} panels x {a.ng} nodes [{time.time()-t0:.0f}s]', flush=True)
    res = {'L': a.L, 'T': a.T, 'K': K, 'w': a.w, 'ng': a.ng, 'b': a.b, 'prec': a.prec, 'comb': [n for n, _, _ in terms],
           'lift': a.lift, **info, 'A': str(AL), 'beta': str(beta), 'M_S': str(MS), 'M0': str(M0)}
    for par in G['pars']:
        t0 = time.time()
        degs = [par + 2 * k for k in range(K)]; dmax = degs[-1]
        Mint = MS / arb.pi() * 2 * L * (2 * dmax + 1) * (2 * L * b).exp()
        Epan = w / 2 * arb(64) / 15 * Mint / ((arb(rho) ** 2 - 1) * arb(rho) ** (2 * (a.ng - 1)))
        Eq = Epan * npan
        p = [(2 * L * (2 * d + 1)).sqrt() * sph_i(d, L / 2) for d in degs]
        sgn = 1 if par == 0 else -1
        A = [[Csum[par][i][j] + sgn * 2 * p[i] * p[j] + (beta if i == j else 0) + arb(0, Eq.upper()) for j in range(K)] for i in range(K)]
        rE = K * max(float(A[i][j].rad()) for i in range(K) for j in range(K))
        Atil = [[arb(A[i][j].mid()) for j in range(K)] for i in range(K)]
        lam_est = lam_estimate(Atil, K)
        print(f' parity {par}: lam_est={lam_est} rE={rE:.3e} quad/entry={float(Eq):.3e} [{time.time()-t0:.0f}s]', flush=True)
        if a.dump_atil:
            import json as _json
            _json.dump([[Atil[i][j].str(radius=False) for j in range(K)] for i in range(K)],
                       open(f"{a.dump_atil}_{'even' if par==0 else 'odd'}.json", 'w'))
            print(f' dumped Atil -> {a.dump_atil}_{"even" if par==0 else "odd"}.json', flush=True)
        lam0 = None
        if lam_est is not None and lam_est > 0 and not a.estimate_only:
            sp = a.schur_prec if a.schur_prec and a.schur_prec > 0 else a.prec
            bslist = [int(x) for x in a.schur_bs.split(',') if x.strip()]
            # rebuild midpoints at schur precision (point balls)
            oldp = ctx.prec; ctx.prec = sp
            Atil_s = [[arb(Atil[i][j].str(radius=False)) for j in range(K)] for i in range(K)]
            print(f' schur_prec={sp} bs_list={bslist}', flush=True)
            for bs in bslist:
                for fac in (0.999, 0.99, 0.95, 0.8, 0.5):
                    c = float('%.4g' % (lam_est * fac)); l0 = arb(repr(c))
                    ok = pd_schur([[Atil_s[i][j] - (l0 if i == j else 0) for j in range(K)] for i in range(K)], bs)
                    print(f'  try bs={bs} fac={fac} l0={c:.6e} -> {ok}', flush=True)
                    if ok:
                        lam0 = l0; break
                if lam0 is not None: break
            ctx.prec = oldp
        # tails (identical to s27)
        dK = par + 2 * K; lt = math.log(Lf * Tf)
        def la(d): return 0.5 * math.log(2 * Lf * (2 * d + 1)) + d * lt - lgdf(d)
        sumD = 0.0; d = dK
        while True:
            term = sexp(2 * la(d))
            r = (Lf * Tf) ** 4 * (2 * d + 5) / ((2 * d + 1) * (2 * d + 3) ** 2 * (2 * d + 5) ** 2)
            if r < 0.5: sumD += 2 * term; break
            sumD += term; d += 2
        sumA = sum(2 * Lf * (2 * dd + 1) for dd in degs); M0f = float(M0)
        lpf = lambda d: 0.5 * math.log(2 * Lf * (2 * d + 1)) + d * math.log(Lf / 2) - lgdf(d) + (Lf / 2) ** 2 / 2
        pD2 = sum(sexp(2 * lpf(dd)) for dd in range(dK, dK + 400, 2)) * 1.0000001
        pA = math.sqrt(sum(float(abs(x).upper()) ** 2 for x in p))
        eps_B = (Tf / math.pi) * M0f * math.sqrt(sumA * sumD) + 2 * pA * math.sqrt(pD2)
        eps_CDD = (Tf / math.pi) * M0f * sumD
        lamD = float(beta.lower()) - eps_CDD - (2 * pD2 if par == 1 else 0)
        out = {'parity': 'even' if par == 0 else 'odd', 'degs': [degs[0], dmax], 'nodes': npan * a.ng, 'quad_err_entry': float(Eq.upper()),
               'lam_est': lam_est, 'normE_bound': rE, 'lam0_schur': (lam0.str(radius=False) if lam0 is not None else None),
               'first_omitted_degree': dK, 'eps_B': eps_B, 'eps_CDD': eps_CDD, 'lamD_lower': lamD}
        if lam0 is not None: out['certified_lower'] = min(float(lam0.lower()) - rE, lamD) - eps_B
        print(json.dumps(out), f'[{time.time()-t0:.0f}s]', flush=True)
        res['even' if par == 0 else 'odd'] = out
    res['wall_s'] = time.time() - T0
    if a.out: json.dump(res, open(a.out, 'w'), indent=1)

if __name__ == '__main__':
    main()
