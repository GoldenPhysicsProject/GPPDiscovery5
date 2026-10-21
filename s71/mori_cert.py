"""Grok S67 CERTIFIED-ASSEMBLY variant of s66/mori_route.py (ball radii kept end to end, exact rational L,s,m,Omega,
delta = 2^-delta_log2 <= d(c) checked, T dumped as midpoint strings + rigorous radius bound for s36_cholres/s37_ldlt_exact).
Original docstring follows.
Grok S66: Mori/Liu bounded-comparison finite condition (Mori VII = [9] (3.4); Mori VIII (2.2)), our own
implementation (does not use Mori's code).  For one parity p:
  A = E*ME, B = E*M^2E with M = gamma I - C_p  (x-space, Gauss-Legendre on cells cut at the shift breakpoints,
      exact for the polynomial integrands);
  J = E*(K+U)E: K = 2|cosh x/2><cosh x/2| - 2|sinh x/2><sinh x/2| + (1/2pi) int_{|t|<Omega} (a(t)-gamma)|e_t><e_t|,
      band part (1/pi) int_0^Omega (a-gamma) b_i b_j dt, b_j = (-1)^((j-p)/2) sqrt(2L(2j+1)) j_j(Lt) (panels, Gauss);
      U = u_p |E*h_p><E*h_p|, h_p = (I - B_band - delta) v_p  (p = 0: v0 const, p = 1: v1 linear);
  X = (B - mA)^{-1}(A - mI) (Loewner minimiser of W),  W(X) = X + X^T - X^T A X + m^{-1}(I - AX - X^T A + X^T B X),
  T = W(X)^{-1} + J - alpha I - beta (B - A^2)   (alpha, beta from mori_consts.py, Lemma 3.2 with rho = 3).
Output: lambda_min(T) estimate (midpoint arithmetic, inverse power iteration): NUMERICAL, not a certificate.
Also lambda_min(A + J) (the simple compression, Mori VII Remark 4.1) for comparison."""
import sys, os, json, time, argparse, math
import multiprocessing as mp
from flint import arb, acb, arb_mat, ctx
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's27')); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's66'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's53'))
from s27_onestroke_arb import sph_j, sph_i
from mori_consts import consts
G = {}
BOOST = 1.75
def sph_j_all(x, D):
    # S67-resumed fix: Arb's ball propagation through bessel_j at large x is ~e^x-crude (radius 1e+175 at x~760, 1e-31 at x~284),
    # so evaluate at the EXACT midpoint and add the input radius times sup|j_d'| <= 1 (|j_n| <= 1 on R and
    # j_n' = (n j_{n-1} - (n+1) j_{n+1})/(2n+1)).  The downward recurrence runs at boosted precision on the exact point.
    rx = arb(x.rad()); x = arb(x.mid())
    old = ctx.prec; ctx.prec = old + int(BOOST * float(x.mid())) + 64
    try:
        v = [None] * (D + 2); v[D + 1] = sph_j(D + 1, x); v[D] = sph_j(D, x)
        for d in range(D, 0, -1): v[d - 1] = (2 * d + 1) / x * v[d] - v[d + 1]
    finally: ctx.prec = old
    return [+t + arb(0, rx) for t in v]

def parseL(s):
    if s.startswith('log'):   # 'log15/2' -> log(15)/2
        a, b = s[3:].split('/'); return arb(int(a)).log() / int(b)
    p, _, q = s.partition('/'); return arb(int(p)) / (int(q) if q else 1)   # S67: exact rational

def mangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            k = n
            while k % p == 0: k //= p
            return math.log(p) if k == 1 else None
    return None

def legv(y, N):   # orthonormal v_d at x = L y, d < N
    # S67-resumed fix: the forward Legendre recurrence in ball arithmetic wraps radii by (|y|+sqrt(1+y^2)) per degree
    # (~1.27 bits/degree at |y|=1: 2^-512 -> inf by N ~ 400).  Run it on the EXACT midpoint at boosted precision
    # (N log2(|y|+sqrt(1+y^2)) + 64 extra bits) and add rad(y) * d(d+1)/2 (Markov |P_d'| <= d(d+1)/2 on [-1,1]).
    assert bool(abs(y) < 1)
    ry = arb(y.rad()); y0 = arb(y.mid()); fy = abs(float(y0))
    old = ctx.prec; ctx.prec = old + int(N * math.log2(fy + math.sqrt(1 + fy * fy))) + 64
    try:
        P = [arb(1), y0]
        for d in range(1, N - 1): P.append(((2 * d + 1) * y0 * P[d] - d * P[d - 1]) / (d + 1))
    finally: ctx.prec = old
    return [(+P[d] + arb(0, ry * (d * (d + 1) // 2))) * G['nrm'][d] for d in range(N)]

def ser(x): return x.str(n=int(ctx.prec * 0.302) + 5, radius=True)   # S67: ball string (contains x)
def sermid(x): return x.mid().str(n=int(ctx.prec * 0.302) + 5, radius=False)

def work_cell(job):
    ctx.prec = G['prec']; ic, ks = job; x0, x1 = G['bpa'][ic], G['bpa'][ic + 1]
    L, N, p = G['L'], G['N'], G['p']; degs = list(range(p, N, 2)); Np = len(degs)
    h = (x1 - x0) / 2; xc = (x0 + x1) / 2
    V = arb_mat(Np, len(ks)); F = arb_mat(Np, len(ks)); FW = arb_mat(len(ks), Np)
    for c, k in enumerate(ks):
        r, w = G['gl'][k]; x = xc + h * r; wk = 2 * h * w     # factor 2: integrate over [0,L] only (even integrand)
        v = legv(x / L, N)
        f = [arb(0)] * N
        for (a, ca) in G['shifts']:
            y = x - a
            if (y > -L) and (y < L):
                va = legv(y / L, N)
                f = [f[d] + ca * va[d] for d in range(N)]
        for i, d in enumerate(degs):
            V[i, c] = v[d]; F[i, c] = f[d]; FW[c, i] = f[d] * wk
    Cp = V * FW; Cp2 = F * FW
    return [[ser(Cp[i, j]) for j in range(Np)] for i in range(Np)], [[ser(Cp2[i, j]) for j in range(Np)] for i in range(Np)]

def work_panel(job):
    ctx.prec = G['prec']; ip, ng = job; t0, t1 = G['edges'][ip], G['edges'][ip + 1]
    L, N, p, g = G['L'], G['N'], G['p'], G['gamma']; degs = list(range(p, N, 2)); Np = len(degs)
    h = (t1 - t0) / 2; tc = (t0 + t1) / 2
    gl = [arb.legendre_p_root(ng, k, weight=True) for k in range(ng)]
    Bm = arb_mat(Np, ng); BW = arb_mat(ng, Np); bb0 = [arb(0)] * Np
    for k, (r, w) in enumerate(gl):
        t = tc + h * r; wk = h * w / arb.pi()
        js = sph_j_all(L * t, N)
        S = acb(arb(1) / 4, t / 2).digamma().real - G['logpi'] - g
        b = [((-1) ** ((d - p) // 2)) * (2 * L * (2 * d + 1)).sqrt() * js[d] for d in degs]
        for i in range(Np):
            Bm[i, k] = b[i]; BW[k, i] = b[i] * wk * S; bb0[i] += b[i] * b[0] * wk
    Ja = Bm * BW
    return [[ser(Ja[i, j]) for j in range(Np)] for i in range(Np)], [ser(x) for x in bb0]

def lam_near0(Mrows, iters=40):
    A = arb_mat(Mrows); n = A.nrows(); Ai = A.inv().mid()
    v = arb_mat([[arb(1) / (i + 1)] for i in range(n)]); lam = None
    for _ in range(iters):
        u = Ai * v
        nu = sum((u[i, 0] ** 2 for i in range(n)), arb(0)).sqrt()
        rq = sum((u[i, 0] * v[i, 0] for i in range(n)), arb(0)) / sum((v[i, 0] ** 2 for i in range(n)), arb(0))
        lam = 1 / rq; v = arb_mat([[(u[i, 0] / nu).mid()] for i in range(n)])
    return lam, v

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', default='log12/2'); ap.add_argument('--N', type=int, default=832)
    ap.add_argument('--Omega', default='416'); ap.add_argument('--s', default='57/16'); ap.add_argument('--m', default='1/2')
    ap.add_argument('--u0', default='207'); ap.add_argument('--u1', default='68')
    ap.add_argument('--parity', type=int, default=0); ap.add_argument('--prec', type=int, default=448)
    ap.add_argument('--workers', type=int, default=8); ap.add_argument('--chunk', type=int, default=64)
    ap.add_argument('--pw', default='16'); ap.add_argument('--png', type=int, default=96)
    ap.add_argument('--delta_log2', type=int, default=0); ap.add_argument('--out', default=''); ap.add_argument('--dump', default='')
    a = ap.parse_args()
    ctx.prec = a.prec; t_start = time.time()
    L = parseL(a.L); N = a.N; p = a.parity; Om = parseL(a.Omega); s = parseL(a.s); m = parseL(a.m); g = s + m
    G.update(prec=a.prec, L=L, N=N, p=p, gamma=g, logpi=arb.pi().log())
    G['nrm'] = [((2 * d + 1) / (2 * L)).sqrt() for d in range(N)]
    e2L = (2 * L).exp(); sh = []
    n = 2
    while arb(n) < e2L:
        lam = mangoldt(n)
        if lam is not None:
            c = arb(n).log() / arb(n).sqrt() if lam == math.log(n) else arb([q for q in range(2, n + 1) if n % q == 0][0]).log() / arb(n).sqrt()
            sh += [(arb(n).log(), c), (-arb(n).log(), c)]
        n += 1
    G['shifts'] = sh
    cst, ok = consts(L, Om, s, m, int(a.u0), int(a.u1), N)
    print('consts ok', ok, {k: v.str(6) for k, v in cst.items()}, flush=True)
    assert a.delta_log2 > 0 and bool(arb(a.delta_log2) > cst['mlog2_dc']), 'need 0 < delta = 2^-k <= d(c) (S69 fix: checked before assembly)'
    # cells on [0, L]
    #bps = sorted(set([0.0, float(L.mid())] + [abs(float(ell.mid())) - float(L.mid()) if float(ell.mid()) > 0 else None for (ell, _) in sh if float(ell.mid()) > 0]) - {None})
    cand = []
    for (ell, _) in sh:
        if float(ell.mid()) > 0:
            z = ell - L if float((ell - L).mid()) > 0 else L - ell
            if 0 < float(z.mid()) < float(L.mid()) and all(abs(float((z - y).mid())) > 1e-30 for y in cand): cand.append(z)
    bpa = [arb(0)] + sorted(cand, key=lambda z: float(z.mid())) + [L]
    G['gl'] = [arb.legendre_p_root(N, k, weight=True) for k in range(N)]
    jobs = []
    for i in range(len(bpa) - 1):
        for k0 in range(0, N, a.chunk): jobs.append((i, list(range(k0, min(N, k0 + a.chunk)))))
    # panels on [0, Omega]: geometric 0,1/2,1,2,..,pw then uniform width pw
    edges = [arb(0), arb(1) / 2, arb(1)]; pw = arb(a.pw)
    while edges[-1] < pw: edges.append(edges[-1] * 2)
    while edges[-1] < Om: edges.append(edges[-1] + pw)
    edges[-1] = Om
    G['bpa'] = bpa; G['edges'] = edges
    pjobs = [(i, a.png) for i in range(len(edges) - 1)]
    print('cells', len(bpa) - 1, 'cell jobs', len(jobs), 'panels', len(pjobs), flush=True)
    Np = len(range(p, N, 2))
    CpA = arb_mat(Np, Np); Cp2A = arb_mat(Np, Np); JaA = arb_mat(Np, Np); bb0 = [arb(0)] * Np
    with mp.get_context('fork').Pool(a.workers) as pool:   # streaming reduction (S66 fix: OOM when holding all partials)
        for (j1, b1) in pool.imap_unordered(work_panel, pjobs):
            JaA += arb_mat([[arb(x) for x in row] for row in j1]); bb0 = [bb0[i] + arb(b1[i]) for i in range(Np)]
        print('panels done %.0fs' % (time.time() - t_start), flush=True)
        for (c1, c2) in pool.imap_unordered(work_cell, jobs):
            CpA += arb_mat([[arb(x) for x in row] for row in c1]); Cp2A += arb_mat([[arb(x) for x in row] for row in c2])
        print('cells done %.0fs' % (time.time() - t_start), flush=True)
    Cp = CpA.tolist(); Cp2 = Cp2A.tolist(); Ja = JaA.tolist()
    degs = list(range(p, N, 2))
    if p == 0: pv = [(2 * L * (2 * d + 1)).sqrt() * sph_i(d, L / 2) for d in degs]; psign = 2
    else: pv = [(2 * L * (2 * d + 1)).sqrt() * sph_i(d, L / 2) for d in degs]; psign = -2
    assert a.delta_log2 > 0 and bool(arb(a.delta_log2) > cst['mlog2_dc']), 'need 0 < delta = 2^-k <= d(c)'
    delta = arb(2) ** (-a.delta_log2)
    hvec = [((1 - delta) if i == 0 else arb(0)) - bb0[i] for i in range(Np)]
    uu = arb(int(a.u0)) if p == 0 else arb(int(a.u1))
    I = arb_mat([[arb(1) if i == j else arb(0) for j in range(Np)] for i in range(Np)])
    CpM = arb_mat(Cp); CpM = (CpM + CpM.transpose()) / 2
    Cp2M = arb_mat(Cp2); Cp2M = (Cp2M + Cp2M.transpose()) / 2
    A = g * I - CpM; B = g * g * I - 2 * g * CpM + Cp2M
    J = arb_mat([[Ja[i][j] + psign * pv[i] * pv[j] + uu * hvec[i] * hvec[j] for j in range(Np)] for i in range(Np)])
    J = (J + J.transpose()) / 2
    def rmax(M): return max(float(M[i, j].rad()) for i in range(Np) for j in range(Np))
    if os.environ.get('MORI_DIAG'):
        print('DIAG rad Cp %.3e Cp2 %.3e Ja %.3e pv %.3e hvec %.3e bb0 %.3e' % (rmax(CpM), rmax(Cp2M), max(float(arb(x).rad()) for r in Ja for x in r), max(float(x.rad()) for x in pv), max(float(x.rad()) for x in hvec), max(float(x.rad()) for x in bb0)), flush=True)
        bad = [(i, j) for i in range(Np) for j in range(Np) if not arb(Ja[i][j]).rad().is_finite()][:5]; print('DIAG nonfinite Ja', bad, flush=True)
        badc = [(i, j) for i in range(Np) for j in range(Np) if not CpM[i, j].rad().is_finite()][:5]; print('DIAG nonfinite Cp', badc, flush=True)
    print('max radii A %.3e B %.3e J %.3e' % (rmax(A), rmax(B), rmax(J)), flush=True)
    # S67: X is ANY matrix (W(X) >= G holds for every X): take the point matrix from midpoints, approx solve
    X = (B - m * A).mid().solve((A - m * I).mid(), algorithm='approx').mid(); XT = X.transpose()
    W = X + XT - XT * A * X + (I - A * X - XT * A + XT * B * X) / m
    W = (W + W.transpose()) / 2
    print('max radius W %.3e' % rmax(W), flush=True)
    # S67: rigorous inverse by residual: Y point approx inverse, R = I - Y W (balls); if ||R|| < 1 then
    # ||W^-1 - Y||_2 <= ||R||_2 ||Y||_2 / (1 - ||R||_2), and ||W^-1 - sym(Y)||_2 <= the same (W^-1 symmetric).
    Y = W.mid().solve(I, algorithm='approx').mid(); Y = ((Y + Y.transpose()) / 2).mid()
    R = I - Y * W
    def n2(M):  # ||M||_2 <= sqrt(||M||_1 ||M||_inf), upper bound from |entries|
        ri = max(sum((abs(M[i, j]).upper() for j in range(Np)), arb(0)) for i in range(Np))
        ci = max(sum((abs(M[i, j]).upper() for i in range(Np)), arb(0)) for j in range(Np))
        return arb((ri * ci).sqrt().upper())
    nR = n2(R); nY = n2(Y); assert bool(nR < 1)
    dinv = arb((nR * nY / (1 - nR)).upper())
    G['dinv'] = dinv; G['nR'] = nR
    print('||R||_2 <= %s, ||Y||_2 <= %s, ||W^-1 - Y||_2 <= %s' % (nR.str(5), nY.str(5), dinv.str(5)), flush=True)
    T = Y + J - cst['alpha'] * I - cst['beta'] * (B - A * A)
    T = (T + T.transpose()) / 2
    K = Np
    radT = max(sum((arb(T[i, j].rad()) for j in range(K)), arb(0)) for i in range(K))   # >= ||T - mid T||_2
    radW = max(sum((arb(W[i, j].rad()) for j in range(K)), arb(0)) for i in range(K))
    rows = [[sermid(T[i, j]) for j in range(K)] for i in range(K)]
    import hashlib
    raw = json.dumps(rows).encode(); open(a.dump, 'wb').write(raw); sha = hashlib.sha256(raw).hexdigest()
    res = dict(L=a.L, N=N, parity=p, Omega=a.Omega, s=a.s, m=a.m, u0=a.u0, u1=a.u1, prec=a.prec, pw=a.pw, png=a.png,
               delta_log2=a.delta_log2, dump=a.dump, sha256=sha, radT_rowsum=float(radT.upper()), radW_rowsum=float(radW.upper()), inv_err=float(G['dinv'].upper()), normR=float(G['nR'].upper()),
               T00=T[0, 0].str(10), consts={k: v.str(8) for k, v in cst.items()}, consts_ok=ok, seconds=round(time.time() - t_start))
    print(json.dumps(res), flush=True)
    if a.out:
        with open(a.out, 'a') as fh: fh.write(json.dumps(res) + '\n')

if __name__ == '__main__':
    main()
