"""Grok S66: Mori/Liu bounded-comparison finite condition (Mori VII = [9] (3.4); Mori VIII (2.2)), our own
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
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's27'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's53'))
from s27_onestroke_arb import sph_j, sph_i
from mori_consts import consts
G = {}
BOOST = 1.75
def sph_j_all(x, D):
    old = ctx.prec; ctx.prec = old + int(BOOST * float(x.mid())) + 64
    try:
        v = [None] * (D + 2); v[D + 1] = sph_j(D + 1, x); v[D] = sph_j(D, x)
        for d in range(D, 0, -1): v[d - 1] = (2 * d + 1) / x * v[d] - v[d + 1]
    finally: ctx.prec = old
    return [+t for t in v]

def parseL(s):
    if s.startswith('log'):   # 'log15/2' -> log(15)/2
        a, b = s[3:].split('/'); return arb(int(a)).log() / int(b)
    return arb(s)

def mangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            k = n
            while k % p == 0: k //= p
            return math.log(p) if k == 1 else None
    return None

def legv(y, N):   # orthonormal v_d at x = L y, d < N
    P = [arb(1), y]
    for d in range(1, N - 1): P.append(((2 * d + 1) * y * P[d] - d * P[d - 1]) / (d + 1))
    return [P[d] * G['nrm'][d] for d in range(N)]

def ser(x): return x.mid().str(n=int(ctx.prec * 0.302) + 5, radius=False)

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
    L = parseL(a.L); N = a.N; p = a.parity; Om = arb(a.Omega); s = arb(a.s); m = arb(a.m); g = s + m
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
    delta = arb(2) ** (-a.delta_log2) if a.delta_log2 else arb(0)
    hvec = [((1 - delta) if i == 0 else arb(0)) - bb0[i] for i in range(Np)]
    uu = arb(int(a.u0)) if p == 0 else arb(int(a.u1))
    I = arb_mat([[arb(1) if i == j else arb(0) for j in range(Np)] for i in range(Np)])
    CpM = arb_mat(Cp); CpM = (CpM + CpM.transpose()) / 2
    Cp2M = arb_mat(Cp2); Cp2M = (Cp2M + Cp2M.transpose()) / 2
    A = g * I - CpM; B = g * g * I - 2 * g * CpM + Cp2M
    J = arb_mat([[Ja[i][j] + psign * pv[i] * pv[j] + uu * hvec[i] * hvec[j] for j in range(Np)] for i in range(Np)])
    J = (J + J.transpose()) / 2
    X = (B - m * A).solve(A - m * I).mid(); XT = X.transpose()   # Loewner-optimal X (see mori_post.py)
    W = X + XT - XT * A * X + (I - A * X - XT * A + XT * B * X) / m
    W = (W + W.transpose()) / 2
    T = W.inv() + J - cst['alpha'] * I - cst['beta'] * (B - A * A)
    T = (T + T.transpose()) / 2
    if a.dump:
        import pickle
        with open(a.dump, 'wb') as fh:
            pickle.dump({k: [[ser(M[i, j]) for j in range(Np)] for i in range(Np)] for k, M in (('T', T), ('A', A), ('B', B), ('J', J))}, fh)
        print('dumped', a.dump, flush=True)
    def safe(M):
        try: return lam_near0(M.mid().tolist())
        except ZeroDivisionError: return arb(0), None
    lamT, vT = safe(T)
    lamAJ, _ = safe(A + J)
    lamDM = arb(0)
    # degree profile of the bottom vector of T
    vv = [float(abs(vT[i, 0]).mid()) for i in range(Np)]
    tot = sum(x * x for x in vv); acc = 0; d99 = None
    for i, x in enumerate(vv):
        acc += x * x
        if d99 is None and acc >= 0.99 * tot: d99 = degs[i]
    res = dict(L=a.L, N=N, parity=p, Omega=a.Omega, s=a.s, m=a.m, u0=a.u0, u1=a.u1, prec=a.prec, pw=a.pw, png=a.png,
               lam_T=lamT.str(12), lam_A_plus_J=lamAJ.str(12), lam_near0_DM=lamDM.str(6), deg99_bottom=d99,
               consts={k: v.str(8) for k, v in cst.items()}, consts_ok=ok, seconds=round(time.time() - t_start))
    print(json.dumps(res), flush=True)
    if a.out:
        with open(a.out, 'a') as fh: fh.write(json.dumps(res) + '\n')

if __name__ == '__main__':
    main()
