"""Grok S27: independent Arb (ball-arithmetic) implementation of the one-stroke
frequency-split reduction of X. Zhu, arXiv:2608.24827 (Thm 1.1), for zeta on the
window supp f in [-L, L].

R(f) = 2|<cosh(x/2),f>|^2 - 2|<sinh(x/2),f>|^2 + (1/pi) int_0^T (Psi - beta)|F|^2 dt + beta ||f||^2
with Psi(t) = Re psi(1/4+it/2) - log pi - sum_{log n < 2L} 2 Lambda(n)/sqrt(n) cos(t log n),
beta = log(T/2pi) - 1/T - A_L.  Q(f) >= R(f) by the envelope Re psi(1/4+it/2) - log pi >= log(t/2pi) - 1/t (t >= 15/4).

Rigour: Gauss-Legendre panels with Arb-enclosed nodes/weights; per-panel Gauss error via
Trefethen's Bernstein-ellipse bound (64/15) M rho^{-2n}/(rho^2-1) with M from Arb evaluation of the
analytically continued symbol on complex boxes covering the strip |Im t| <= b; Cholesky of
A - lam0 I in ball arithmetic; explicit super-exponential tail (block coupling) bounds.
Status: computer-assisted (Arb), NOT Lean.  RH OPEN.
"""
import sys, math, json, time, argparse
from flint import arb, acb, arb_mat, ctx

def mangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return math.log(p) if m == 1 else 0.0, p
    return 0.0, None

def comb_terms(L):
    terms = []
    n = 2
    while math.log(n) < 2 * L:
        lam, p = mangoldt(n)
        if lam > 0:
            terms.append((n, p))
        n += 1
    # exactness check of the boundary
    for (nn, p) in terms:
        assert arb(nn).log() < 2 * arb(L)
    assert arb(n).log() > 2 * arb(L)
    return [(n, 2 * arb(p).log() / arb(n).sqrt(), arb(n).log()) for (n, p) in terms]

def sph_j(d, x):
    return (arb.pi() / (2 * x)).sqrt() * x.bessel_j(arb(d) + arb(1) / 2)

def sph_j_all(x, D):
    """j_0..j_{D+1}(x), x>0: top two orders by Arb Bessel at boosted precision (0F1 cancellation ~ e^x),
    then exact backward recurrence j_{d-1} = (2d+1)/x j_d - j_{d+1} in ball arithmetic."""
    old = ctx.prec
    ctx.prec = old + int(1.5 * float(x.mid())) + 64
    try:
        v = [None] * (D + 2)
        v[D + 1] = sph_j(D + 1, x); v[D] = sph_j(D, x)
        for d in range(D, 0, -1):
            v[d - 1] = (2 * d + 1) / x * v[d] - v[d + 1]
    finally:
        ctx.prec = old
    return [+t for t in v]

def sph_i(d, x):
    return (arb.pi() / (2 * x)).sqrt() * x.bessel_i(arb(d) + arb(1) / 2)

def lgdf(d):  # log (2d+1)!!
    return math.lgamma(2 * d + 2) - d * math.log(2) - math.lgamma(d + 1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', default='1/2')
    ap.add_argument('--T', default='35')
    ap.add_argument('--K', type=int, default=30)     # modes per sector
    ap.add_argument('--w', default='1/2')            # panel width
    ap.add_argument('--ng', type=int, default=24)    # Gauss nodes per panel
    ap.add_argument('--b', default='2/5')            # strip half-height
    ap.add_argument('--prec', type=int, default=192)
    ap.add_argument('--parity', default='both')
    ap.add_argument('--out', default=None)
    ap.add_argument('--bis', type=int, default=24)
    ap.add_argument('--lift', default=None)
    ap.add_argument('--kgrid', type=int, default=200000)
    ap.add_argument('--lo', type=float, default=-40.0)
    a = ap.parse_args()
    ctx.prec = a.prec
    L = arb(a.L); T = arb(a.T); w = arb(a.w); b = arb(a.b)
    Lf, Tf, wf, bf = float(L.mid()), float(T.mid()), float(w.mid()), float(b.mid())
    terms = comb_terms(Lf)
    AL = sum((c for (_, c, _) in terms), arb(0))
    lift_terms = []
    kappas = {}
    if a.lift:
        lj = json.load(open(a.lift))
        AL = arb(0)
        for ps, d in lj['primes'].items():
            p = int(ps); lp = arb(p).log()
            ks = d['ks']; cs = [2 * lp / arb(p) ** (arb(k) / 2) for k in ks]
            inv = d['inv']; av = [arb('%.12e' % x) for x in d['a']]
            for k in inv:
                assert arb(k) * lp >= 2 * L, (p, k)
            # rigorous upper bound of kappa_p = max_theta (B_p - h_p) on [0, pi]
            ng = a.kgrid; hw = arb.pi() / (2 * ng); kap = None
            for i in range(ng):
                th = arb(((2 * i + 1) * hw).mid(), (hw * 1.0001).mid())
                v = sum((c * (k * th).cos() for k, c in zip(ks, cs)), arb(0)) - sum((x * (k * th).cos() for k, x in zip(inv, av)), arb(0))
                u = v.upper()
                kap = u if kap is None or u > kap else kap
            kap = arb(kap)
            kappas[p] = float(kap.mid())
            AL += kap
            lift_terms += [(x, k * lp) for k, x in zip(inv, av)]
        print('lift kappas (rigorous upper):', kappas, 'A_lift', AL.mid(), flush=True)
    beta = (T / (2 * arb.pi())).log() - 1 / T - AL
    assert T > arb(15) / 4 and beta > 0
    print(f'L={a.L} T={a.T} comb n={[n for n,_,_ in terms]} A_L={AL.mid()} beta={beta.mid()}', flush=True)
    npan = int(round(Tf / wf)); assert abs(npan * wf - Tf) < 1e-12
    gl = [arb.legendre_p_root(a.ng, k, weight=True) for k in range(a.ng)]
    # Bernstein rho for the strip: (w/2)(rho - 1/rho)/2 = b
    s = 4 * bf / wf; rho = (s + math.sqrt(s * s + 4)) / 2 * (1 - 1e-9)
    ovh = (wf / 2) * ((rho + 1 / rho) / 2 - 1)
    # --- M_S: sup |Psi_analytic - beta| over [-ovh, T+ovh] x [-b, b] via Arb boxes
    t0 = time.time()
    # Max-modulus principle: the continued symbol is analytic on the closed rectangle
    # [-ovh, T+ovh] x [-b, b] (digamma poles sit at |Im z| >= 1/2), so its sup is attained on the boundary.
    hb = 0.002
    MS = arb(0)
    logpi = arb.pi().log()
    def symb(z):
        v = ((acb(arb(1)/4) + acb(0, 1) * z / 2).digamma() + (acb(arb(1)/4) - acb(0, 1) * z / 2).digamma()) / 2 - logpi - beta
        for (_, c, lg) in terms:
            v -= c * (z * lg).cos()
        for (x, fr) in lift_terms:
            v += x * (z * fr).cos()
        return abs(v).upper()
    x0 = -ovh; x1 = Tf + ovh
    nx = int(math.ceil((x1 - x0) / (2 * hb)))
    for ix in range(nx):
        xc = arb(x0) + (2 * ix + 1) * arb(hb)
        for yy in (b, -b):
            m = symb(acb(arb(xc, hb * 1.001), yy))
            if m > MS: MS = m
    ny = int(math.ceil(2 * bf / (2 * hb)))
    for iy in range(ny):
        yc = -b + (2 * iy + 1) * arb(hb)
        for xx in (x0, x1):
            m = symb(acb(arb(xx), arb(yc, hb * 1.001)))
            if m > MS: MS = m
    MS = arb(MS.upper())
    print(f'M_S (strip sup) <= {MS.mid()}  [{time.time()-t0:.1f}s] rho={rho:.4f} ovh={ovh:.3f}', flush=True)
    res = {'L': a.L, 'T': a.T, 'K': a.K, 'w': a.w, 'ng': a.ng, 'b': a.b, 'prec': a.prec,
           'comb': [n for n, _, _ in terms], 'lift': a.lift, 'kappas': kappas, 'A_L': str(AL), 'beta': str(beta), 'M_S': str(MS)}
    pars = [0, 1] if a.parity == 'both' else [int(a.parity)]
    for par in pars:
        t0 = time.time()
        degs = [par + 2 * k for k in range(a.K)]
        K = a.K
        # quadrature: accumulate C = sum_nodes wt * v v^T in chunks (memory-light)
        C = arb_mat(K, K); N = 0
        chunk_cols = []; chunk_w = []
        def flush():
            nonlocal C
            n = len(chunk_cols)
            if n == 0: return
            F = arb_mat(K, n); FW = arb_mat(n, K)
            for j in range(n):
                cj = chunk_cols[j]; wj = chunk_w[j]
                for i in range(K):
                    F[i, j] = cj[i]; FW[j, i] = cj[i] * wj
            C = C + F * FW
            chunk_cols.clear(); chunk_w.clear()
        for ip in range(npan):
            lo = ip * w
            for (xk, wk) in gl:
                t = lo + w / 2 * (1 + xk)
                S = acb(arb(1)/4, t / 2).digamma().real - logpi - beta
                for (_, c, lg) in terms:
                    S -= c * (t * lg).cos()
                for (x, fr) in lift_terms:
                    S += x * (t * fr).cos()
                x = L * t
                jall = sph_j_all(x, degs[-1])
                chunk_cols.append([((-1) ** (d // 2)) * (2 * L * (2 * d + 1)).sqrt() * jall[d] for d in degs])
                chunk_w.append(w / 2 * wk * S / arb.pi()); N += 1
            if len(chunk_cols) >= 2048: flush()
        flush()
        # Gauss error bound per entry (uniform over the block)
        dmax = degs[-1]
        Mint = MS / arb.pi() * 2 * L * (2 * dmax + 1) * (2 * L * b).exp()
        Epan = w / 2 * arb(64) / 15 * Mint / ((arb(rho) ** 2 - 1) * arb(rho) ** (2 * a.ng))
        Eq = Epan * npan
        # pole vector
        p = [(2 * L * (2 * d + 1)).sqrt() * sph_i(d, L / 2) for d in degs]
        sgn = 1 if par == 0 else -1
        A = arb_mat(K, K)
        for i in range(K):
            for j in range(K):
                v = C[i, j] + sgn * 2 * p[i] * p[j] + (beta if i == j else 0) + arb(0, Eq.upper())
                A[i, j] = v
        # Certification strategy (as Zhu Sec. 5): Atil = exact midpoint matrix (zero radius);
        # lam_min(A) >= lam_min(Atil) - ||A - Atil||, ||A - Atil|| <= K * max entry radius.
        Atil = [[arb(A[i, j].mid()) for j in range(K)] for i in range(K)]
        rE = K * max(float(A[i, j].rad()) for i in range(K) for j in range(K))
        def chol_ok(lam0):
            M = [[Atil[i][j] - (lam0 if i == j else 0) for j in range(K)] for i in range(K)]
            for k in range(K):
                for j in range(k):
                    M[k][k] -= M[k][j] * M[k][j]
                if not (M[k][k] > 0):
                    return False
                M[k][k] = M[k][k].sqrt()
                for i in range(k + 1, K):
                    acc = M[i][k]
                    for j in range(k):
                        acc -= M[i][j] * M[k][j]
                    M[i][k] = acc / M[k][k]
            return True
        # bisection in log10 for the largest passing shift
        lo, hi = a.lo, 0.0
        if not chol_ok(arb(10) ** int(lo)):
            lam0 = None; lam_est = None
        else:
            for _ in range(a.bis):
                mid = (lo + hi) / 2
                if chol_ok(arb(10) ** arb(mid)): lo = mid
                else: hi = mid
            # report a short decimal shift that passes
            lam_est = 10 ** hi
            c = float('%.3g' % (10 ** lo))
            lam0 = arb(repr(c))
            while not chol_ok(lam0):
                c = float('%.3g' % (c * 0.99)); lam0 = arb(repr(c))
        evs = None
        # tails
        dK = par + 2 * K
        lt = math.log(Lf * Tf)
        def la(d): return 0.5 * math.log(2 * Lf * (2 * d + 1)) + d * lt - lgdf(d)
        sumD = 0.0; d = dK
        while True:
            term = math.exp(2 * la(d))
            r = (Lf * Tf) ** 4 * (2 * d + 5) / ((2 * d + 1) * (2 * d + 3) ** 2 * (2 * d + 5) ** 2)
            if r < 0.5:
                sumD += 2 * term; break
            sumD += term; d += 2
        sumA = sum(2 * Lf * (2 * dd + 1) for dd in degs)
        M0 = float(MS.upper())
        lp = lambda d: 0.5 * math.log(2 * Lf * (2 * d + 1)) + d * math.log(Lf / 2) - lgdf(d) + (Lf / 2) ** 2 / 2
        pD2 = sum(math.exp(2 * lp(dd)) for dd in range(dK, dK + 400, 2)) * 1.0000001
        pA = math.sqrt(sum(float(abs(x).upper()) ** 2 for x in p))
        eps_B = (Tf / math.pi) * M0 * math.sqrt(sumA * sumD) + 2 * pA * math.sqrt(pD2)
        eps_CDD = (Tf / math.pi) * M0 * sumD
        lamD = float(beta.lower()) - eps_CDD - (2 * pD2 if par == 1 else 0)
        out = {'parity': 'even' if par == 0 else 'odd', 'degs': [degs[0], degs[-1]], 'nodes': N,
               'quad_err_entry': float(Eq.upper()), 'lam_fail_shift': lam_est, 'normE_bound': rE,
               'lam0_cholesky': (lam0.str(radius=False) if lam0 is not None else None),
               'first_omitted_degree': dK, 'sum_tail_amp2': sumD, 'eps_B': eps_B, 'eps_CDD': eps_CDD,
               'lamD_lower': lamD, 'max_entry_radius': max(float(A[i, j].rad()) for i in range(K) for j in range(K))}
        if lam0 is not None:
            out['certified_lower'] = min(float(lam0.lower()) - rE, lamD) - eps_B
        print(json.dumps(out), f'[{time.time()-t0:.1f}s]', flush=True)
        res['even' if par == 0 else 'odd'] = out
    if a.out:
        json.dump(res, open(a.out, 'w'), indent=1)

if __name__ == '__main__':
    main()
