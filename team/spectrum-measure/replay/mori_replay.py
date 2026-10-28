"""Independent replay of the Mori-route finite condition (Mori, 'Weil positivity on the window of width log 12',
Zenodo 10.5281/zenodo.23265728 v2, Sections 3-5), written from the paper for the GPP spectrum-measure seat.
It does not import or copy GPPDiscovery5 team code (s66/s67/s71/s37); conventions were read only.

Objects (one parity p, H = L^2(-L,L), v_d = sqrt((2d+1)/(2L)) P_d(x/L), d<N, d = p mod 2):
  C_p = sum_{n<e^{2L}} Lambda(n) n^{-1/2} (S_{log n} + S_{-log n}),  ||C_p|| <= s (input, certified separately: cp_norm_check.py)
  gamma = s+m, M = gamma - C_p, A = E*ME = gamma - <v,C_p v>, B = E*M^2E = gamma^2 - 2 gamma <v,C_p v> + <C_p v, C_p v>
  J = E*(K+U)E,  K = +-2|pole><pole| + (1/pi) int_0^Omega (a(t)-gamma) b_i b_j dt,  b_d(t) = (-1)^((d-p)/2) sqrt(2L(2d+1)) j_d(Lt)
  U = u_p |h><h|, h = (1-delta) e_p - (1/pi) int_0^Omega b_i b_p dt
  T = W(X)^{-1} + J - alpha - beta (B - A^2),  W(X) = X+X'-X'AX + (I-AX-X'A+X'BX)/m
Differences from the team implementation (by design):
  * A,B: Gauss-Legendre on cells of [0,L] (same mathematics), node-vectorised Legendre recurrence; Markov radius bound added
    as a scalar budget, not per entry.
  * J: adaptive panels WITH a rigorous Gauss-Legendre error bound (Trefethen ATAP Thm 19.3 on Bernstein ellipses,
    |g| bounded by Arb digamma on ball covers of the ellipse; |b_i b_j| <= 2L e^{2L|Im z|} by Cauchy-Schwarz).
  * certificate: unit-lower-triangular congruence Z (recursive block LDL^T at mid precision) and a rigorous
    ball product S = Z (T - l I) Z^T; diagonal dominance of S  =>  T - l I > 0.  Control: a certified negative
    diagonal entry of S_c  =>  explicit vector with v'(T - l_c)v < 0.
Checkpoints: every cell / panel group is a separate file; reruns skip finished pieces.
Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, os, json, time, math, argparse, pickle, hashlib
from fractions import Fraction
import numpy as np
from flint import arb, acb, arb_mat, ctx

WP = 448            # working precision (bits) for matrices
S_FIX = 400         # fixed-point scale for checkpoint files

# ---------------------------------------------------------------- parameters
def parse_num(s):
    """'71/50' -> exact rational arb; 'log12/2' -> log(12)/2 (returns arb, and a Fraction Q if L = log(Q)/2)."""
    if s.startswith('log'):
        a, b = s[3:].split('/'); assert b == '2'
        return arb(int(a)).log() / 2, Fraction(int(a))
    p, _, q = s.partition('/')
    return arb(int(p)) / (int(q) if q else 1), None

def mangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            k = n
            while k % p == 0: k //= p
            return p if k == 1 else None
    return None

def setup(a):
    ctx.prec = WP
    L, LQ = parse_num(a.L); Om, _ = parse_num(a.Omega); s, _ = parse_num(a.s); m, _ = parse_num(a.m)
    g = s + m
    e2L = (2 * L).exp()
    shifts = []          # (ell, c, n)
    n = 2
    while True:
        if LQ is not None:
            if n >= LQ: break
        elif not bool(arb(n) < e2L):
            assert bool(arb(n) > e2L), 'n = e^{2L} ambiguous'
            break
        p = mangoldt(n)
        if p is not None:
            shifts.append((arb(n).log(), arb(p).log() / arb(n).sqrt(), n))
        n += 1
    return dict(L=L, LQ=LQ, Om=Om, s=s, m=m, g=g, shifts=shifts, N=a.N, p=a.parity)

def a_of_t(t):   # a(t) = Re psi(1/4 + i t/2) - log pi
    return acb(arb(1) / 4, t / 2).digamma().real - arb.pi().log()

def consts(P, a):
    """Mori [VII] Prop. 3.1 and Lemma 3.2 constants, recomputed from the formulas (rho = 3)."""
    L, Om, s, m, g, N = P['L'], P['Om'], P['s'], P['m'], P['g'], P['N']
    out = {}
    c = Om * L
    C_log = (Om / (2 * arb.pi())).log() - 1 / Om - g         # (3.1), uses Mori[14] Lemma 3.3 (a(t) >= log(t/2pi) - 1/t, t >= 15/4)
    C_mono = a_of_t(Om) - g                                    # alternative: a increasing on t >= 0 (elementary), so a(t)-gamma >= a(Omega)-gamma
    dc = (-2 * (c + 1)).exp() / (4 * (c + 1) * (c + 4))
    mlog2_dc = -dc.log() / arb(2).log()
    E0 = (c + 1) / (arb.pi() * c * c)
    E1 = (3 * c * c + 6 * c + 2) / (arb.pi() * c ** 3)
    u0, u1 = arb(a.u0), arb(a.u1)
    rho = arb(3)
    rN2 = 8 * L * rho ** 2 / (rho - 1) ** 2 * (c * (rho - 1 / rho)).exp() * rho ** (-2 * N)
    rN = rN2.sqrt()
    K0 = (Om / arb.pi()) * arb.max(abs(a_of_t(arb(0)) - g), abs(a_of_t(Om) - g))
    Hp = (2 * L).sqrt() * (L / 2).cosh()
    qN = (2 * L).sqrt() * (L / 2).exp() * (L / 2) ** N / arb.fac_ui(N)
    # e: off-diagonal block. Derived here: ||(I-P)U P|| <= sum u_j ||(I-P)h_j|| ||h_j|| <= sqrt(Omega/pi) r_N sum u_j sqrt(E_j)
    e = (K0 * (2 * L).sqrt() + (Om / arb.pi()).sqrt() * (u0 * E0.sqrt() + u1 * E1.sqrt())) * rN + 4 * Hp * qN
    e_paper_literal = (K0 * (2 * L).sqrt() + (Om / arb.pi()).sqrt() * (u0 * E0 + u1 * E1)) * rN + 4 * Hp * qN
    nn = K0 * rN2 + 4 * qN ** 2
    h = (K0 + (Om / arb.pi()) * (u0 + u1)) * rN2 + 4 * qN ** 2
    alpha = e + 2 * e * e / (m - nn)
    beta = (e + nn) / m ** 2 + 2 * h * h / (m ** 2 * (m - nn))
    delta_log2 = a.delta_log2
    checks = dict(C_pos=bool(C_log > 0), C_mono_ge_C_log=bool(C_mono >= C_log), c_gt_2=bool(c > 2), Om_ge_15_4=bool(Om >= arb(15) / 4),
                  delta_le_dc=bool(arb(delta_log2) > mlog2_dc), u0E0_le_C=bool(u0 * E0 <= C_log), u1E1_le_C=bool(u1 * E1 <= C_log),
                  n_lt_m=bool(nn < m))
    out.update(gamma=g, b=2 * s + m, c=c, C=C_log, C_mono=C_mono, mlog2_dc=mlog2_dc, E0=E0, E1=E1, u0E0=u0 * E0, u1E1=u1 * E1,
               rN=rN, K0=K0, Hp=Hp, qN=qN, e=e, e_paper_literal=e_paper_literal, n=nn, h=h, alpha=alpha, beta=beta,
               final_constant_log2=(C_log.log() / arb(2).log() - delta_log2))
    return out, checks

# ---------------------------------------------------------------- fixed-point checkpoint I/O
def mat_to_fixed(M):
    rows, cols = M.nrows(), M.ncols(); ints = []; rmax = 0.0
    for i in range(rows):
        for j in range(cols):
            x = M[i, j]; man, ex = x.mid().man_exp(); man = int(man); ex = int(ex) + S_FIX
            ints.append(man << ex if ex >= 0 else man >> (-ex))
            r = float(x.rad())
            if r > rmax: rmax = r
    return dict(rows=rows, cols=cols, ints=ints, rad=rmax * (1 + 1e-12) + 2.0 ** (-S_FIX))

def fixed_to_mat(D, extra_rad=0.0):
    old = ctx.prec; ctx.prec = WP + 64
    sc = arb(2) ** (-S_FIX); R = D['rad'] + extra_rad
    it = iter(D['ints'])
    M = arb_mat([[arb(next(it)) * sc for j in range(D['cols'])] for i in range(D['rows'])])
    ctx.prec = old
    M = +M if False else M
    Rm = arb_mat([[arb(0, R)] * D['cols'] for _ in range(D['rows'])])
    return M + Rm

def save(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'wb') as fh: pickle.dump(obj, fh, protocol=4)
    os.replace(tmp, path)

# ---------------------------------------------------------------- Legendre values at many points (vectorised)
def legendre_rows(ys, N, p):
    """ys: list of exact (dyadic) arb points in (-1,1).  Returns dict d -> numpy object array of P_d(ys) (balls),
    for d = p mod 2, d < N, computed by the three-term recurrence in ball arithmetic at boosted precision."""
    ymax = max(abs(float(y)) for y in ys)
    old = ctx.prec
    ctx.prec = WP + int(N * math.log2(ymax + math.sqrt(1 + ymax * ymax))) + 96
    try:
        Y = np.array(ys, dtype=object)
        P0 = np.array([arb(1)] * len(ys), dtype=object); P1 = Y.copy()
        rows = {}
        if p == 0: rows[0] = P0
        if p == 1 and N > 1: rows[1] = P1
        for d in range(1, N - 1):
            P2 = ((2 * d + 1) * Y * P1 - d * P0) / (d + 1)
            P0, P1 = P1, P2
            if (d + 1) % 2 == p: rows[d + 1] = P2
    finally:
        ctx.prec = old
    return {d: np.array([+x for x in v], dtype=object) for d, v in rows.items()}

# ---------------------------------------------------------------- cells: A, B pieces
def breakpoints(P):
    L, LQ = P['L'], P['LQ']; cand = []
    for (ell, c, n) in P['shifts']:
        z = ell - L if bool(ell > L) else L - ell
        assert bool(z > 0) or bool(z < 0) or True
        if not (bool(z > 0) and bool(z < L)): 
            assert bool(z >= L) or bool(z <= 0), 'undecided breakpoint'
            continue
        key = (Fraction(n * n) / LQ if LQ is not None else None)
        if key is not None and key < 1: key = 1 / key
        cand.append((z, key))
    pts = []
    for z, key in sorted(cand, key=lambda t: float(t[0].mid())):
        if pts and ((key is not None and key == pts[-1][1])):
            continue
        if pts and key is None:
            assert not z.overlaps(pts[-1][0]), 'breakpoints overlap'
        pts.append((z, key))
    return [arb(0)] + [z for z, _ in pts] + [L]

GL = {}
def cell_job(ic):
    P, chunk, outdir = GL['P'], GL['chunk'], GL['out']; x0, x1 = GL['bpa'][ic], GL['bpa'][ic + 1]
    path = os.path.join(outdir, 'cells', 'cell_%02d.pkl' % ic)
    if os.path.exists(path): return ic, 'cached'
    ctx.prec = WP; t0 = time.time()
    L, N, p = P['L'], P['N'], P['p']; degs = list(range(p, N, 2)); Np = len(degs)
    nrm = [((2 * d + 1) / (2 * L)).sqrt() for d in degs]
    h = (x1 - x0) / 2; xc = (x0 + x1) / 2
    # indicator pattern of shifted points: constant on the open cell; decide at the centre
    act = []
    for (ell, c, n) in P['shifts']:
        for sg in (1, -1):
            y = xc - sg * ell
            inside = bool(y > -L) and bool(y < L); outside = bool(y < -L) or bool(y > L)
            assert inside or outside, 'centre on a breakpoint'
            if inside: act.append((sg * ell, c))
    K = N   # Gauss nodes: integrand degree <= 2N-2
    CpA = arb_mat(Np, Np); Cp2A = arb_mat(Np, Np)
    ry_max = 0.0
    for k0 in range(0, K, chunk):
        ks = list(range(k0, min(K, k0 + chunk))); nk = len(ks)
        nodes = []
        for k in ks:
            r, w = arb.legendre_p_root(K, k, weight=True)
            nodes.append((xc + h * r, 2 * h * w))       # factor 2: integrate over [-L,L] by evenness
        ys = []; 
        for (x, w) in nodes:
            y = x / L; ys.append(arb(y.mid())); ry_max = max(ry_max, float(y.rad()))
        for (a_, c) in act:
            for (x, w) in nodes:
                y = (x - a_) / L; ys.append(arb(y.mid())); ry_max = max(ry_max, float(y.rad()))
        rows = legendre_rows(ys, N, p)
        V = arb_mat(Np, nk); FW = arb_mat(nk, Np); F = arb_mat(Np, nk)
        cs = [c for (_, c) in act]
        for i, d in enumerate(degs):
            r = rows[d]
            base = r[:nk]
            f = [arb(0)] * nk
            for jj, c in enumerate(cs):
                seg = r[nk * (jj + 1): nk * (jj + 2)]
                for q in range(nk): f[q] = f[q] + c * seg[q]
            for q in range(nk):
                V[i, q] = base[q] * nrm[i]; fv = f[q] * nrm[i]; F[i, q] = fv; FW[q, i] = fv * nodes[q][1]
        CpA += V * FW; Cp2A += F * FW
        del rows
    # scalar budget for node uncertainty: |P_d(y)-P_d(y0)| <= |y-y0| d(d+1)/2 (Markov), sup|v_d| <= sqrt((2d+1)/2L)
    csum = sum(float(c.upper()) for (_, c) in act) * (1 + 1e-12)
    dmax = degs[-1]; vmax = math.sqrt((2 * dmax + 1) / (2 * float(L.lower()))) * (1 + 1e-12)
    epsV = ry_max * dmax * (dmax + 1) / 2 * vmax; epsF = epsV * csum
    Fmax = vmax * csum
    wsum = 2 * float(h.upper()) * 2 * (1 + 1e-9)
    extra_Cp = wsum * (epsV * Fmax + vmax * epsF + epsV * epsF)
    extra_Cp2 = wsum * (2 * epsF * Fmax + epsF * epsF)
    save(path, dict(Cp=mat_to_fixed(CpA), Cp2=mat_to_fixed(Cp2A), extra_Cp=extra_Cp, extra_Cp2=extra_Cp2,
                    n_active_shifts=len(act), seconds=time.time() - t0, ry_max=ry_max))
    return ic, 'done %.0fs' % (time.time() - t0)

# ---------------------------------------------------------------- panels: J band part and the h-vector
def rho_of(w):
    """Bernstein-ellipse parameter of complex point w (acb) relative to [-1,1]: |w + sqrt(w^2-1)| (larger branch)."""
    sq = (w * w - 1).sqrt(); r1 = abs(w + sq); r2 = abs(w - sq)
    return arb.max(r1, r2)

def panel_plan(P, target_log10=-140):
    """Adaptive panels on [0, Omega] with rigorous Gauss error bound per entry <= 10^target / (#panels)."""
    L, Om, g = P['L'], P['Om'], P['g']
    edges = [arb(0), arb(1) / 4, arb(1) / 2, arb(1), arb(2), arb(4), arb(8)]
    while bool(edges[-1] + 16 < Om): edges.append(edges[-1] + 16)
    edges.append(Om)
    plan = []
    for i in range(len(edges) - 1):
        t0, t1 = edges[i], edges[i + 1]; h = (t1 - t0) / 2; c = (t0 + t1) / 2
        # nearest poles of g: z = +-i(2k+1/2)
        rp = None
        for k in range(0, 200):
            for sg in (1, -1):
                z = acb(arb(0), sg * (2 * k + arb(1) / 2)); r = rho_of((z - c) / h)
                rp = r if rp is None else arb.min(rp, r)
        u_pole = float(rp.lower().log()) * 0.85
        best = None
        for n in (32, 48, 64, 80, 96, 128, 160, 192, 256):
            # choose u to minimise 2L h sinh(u) - 2 n u subject to u <= u_pole
            u_opt = math.acosh(max(1.0, n / (float(L) * float(h))))
            u = min(u_opt, u_pole)
            est = (2 * float(L) * float(h) * math.sinh(u) - 2 * n * u) / math.log(10)
            if est < target_log10 - 6:
                best = (n, u); break
        assert best is not None, 'panel %d needs more nodes' % i
        plan.append((i, t0, t1, best[0], best[1]))
    return plan

def ellipse_gmax(P, c, h, rho, M=1024):
    for MM in (M, 2048, 8192, 32768):
        r = _ellipse_gmax(P, c, h, rho, MM)
        if r is not None: return r
    raise AssertionError('pole on ellipse cover')

def psi_ball(w):
    """Rigorous digamma on a complex ball: psi(w) = psi(w+K) - sum_{k<K} 1/(w+k) with Re(w+K) >= 10
    (Arb's acb digamma returns nan on wide balls with Re w < 0 at high precision)."""
    K = max(0, math.ceil(-float(w.real.mid())) + 10)
    s = (w + K).digamma()
    for k in range(K): s -= 1 / (w + k)
    return s

def _ellipse_gmax(P, c, h, rho, M):
    old = ctx.prec; ctx.prec = 128
    try: return __ellipse_gmax(P, c, h, rho, M)
    finally: ctx.prec = old

def __ellipse_gmax(P, c, h, rho, M):
    """Rigorous upper bound of |g(z)|, g = (psi(1/4+iz/2)+psi(1/4-iz/2))/2 - log pi - gamma, on the ellipse boundary
    (hence inside, by the maximum principle; holomorphy checked by pole distances)."""
    g = P['g']; lp = arb.pi().log(); best = arb(0)
    for k in range(M):
        th0 = 2 * arb.pi() * k / M; th1 = 2 * arb.pi() * (k + 1) / M; thc = (th0 + th1) / 2
        z = c + h * (rho * acb(thc.cos(), thc.sin()) + (1 / rho) * acb(thc.cos(), -thc.sin())) / 2
        arcrad = h * (rho + 1 / rho) / 2 * (th1 - th0) / 2      # |dz/dth| <= h(rho+1/rho)/2
        zb = acb(arb(z.real.mid(), arcrad.upper() + z.real.rad()), arb(z.imag.mid(), arcrad.upper() + z.imag.rad()))
        w1 = acb(arb(1) / 4) + acb(0, 1) * zb / 2; w2 = acb(arb(1) / 4) - acb(0, 1) * zb / 2
        val = (psi_ball(w1) + psi_ball(w2)) / 2 - lp - g
        v = abs(val)
        if not v.is_finite(): return None
        best = arb.max(best, arb(v.upper()))
    return best

def sph_j_rows(xs, N, p):
    """j_d(x) for d < N, d = p mod 2, at exact points xs > 0: downward recurrence from Arb-evaluated top values."""
    xmax = max(float(x) for x in xs)
    old = ctx.prec; ctx.prec = WP + int(1.0 * xmax) + 128
    try:
        X = np.array(xs, dtype=object)
        tgt = ctx.prec
        def top1(x, d):
            # Arb's J_nu(x) loses all accuracy near the turning point x ~ nu at modest precision (observed:
            # x = 1420, nu = 2376.5 at 2000 bits gives a ball of width 1e380).  Raise precision until the
            # enclosure has relative radius <= 2^-(tgt+32), then round back.  Still a rigorous ball throughout.
            pr = tgt
            while True:
                o2 = ctx.prec; ctx.prec = pr
                try:
                    v = (arb.pi() / (2 * x)).sqrt() * x.bessel_j(arb(d) + arb(1) / 2)
                finally:
                    ctx.prec = o2
                if v.is_finite() and v.mid() != 0 and bool(v.rad() < abs(v.mid()) * arb(2) ** (-(tgt + 32))):
                    return +v
                pr *= 2
                assert pr <= 64 * tgt + 100000, 'bessel_j precision blow-up'
        def top(d): return np.array([top1(x, d) for x in xs], dtype=object)
        jn1 = top(N); jn = top(N - 1); rows = {}
        if (N - 1) % 2 == p: rows[N - 1] = jn
        Xi = np.array([1 / x for x in xs], dtype=object)
        for d in range(N - 1, 0, -1):
            jm = (2 * d + 1) * Xi * jn - jn1
            jn1, jn = jn, jm
            if (d - 1) % 2 == p: rows[d - 1] = jm
    finally:
        ctx.prec = old
    return {d: np.array([+x for x in v], dtype=object) for d, v in rows.items()}

def panel_job(gid):
    P, items, outdir = GL['P'], GL['groups'][gid], GL['out']
    path = os.path.join(outdir, 'panels', 'pg_%02d.pkl' % gid)
    if os.path.exists(path): return gid, 'cached'
    ctx.prec = WP; t0 = time.time()
    L, N, p, g = P['L'], P['N'], P['p'], P['g']; degs = list(range(p, N, 2)); Np = len(degs)
    lp = arb.pi().log()
    Ja = arb_mat(Np, Np); bb = [arb(0)] * Np
    errJ = 0.0; errb = 0.0; info = []
    for (i, ta, tb, n, u) in items:
        h = (tb - ta) / 2; c = (ta + tb) / 2
        rho = arb(u).exp()
        gmax = ellipse_gmax(P, c, h, rho)
        imz = h * (rho - 1 / rho) / 2
        Mb = 2 * L * (2 * L * imz).exp()                       # |b_i b_j| on the ellipse
        fac = h / arb.pi() * arb(64) / 15 / ((rho * rho - 1) * rho ** (2 * n))
        eJ = fac * gmax * Mb; eb = fac * Mb
        errJ += float(eJ.upper()) * (1 + 1e-9); errb += float(eb.upper()) * (1 + 1e-9)
        nodes = [arb.legendre_p_root(n, k, weight=True) for k in range(n)]
        ts = [c + h * r for (r, w) in nodes]
        xs = [arb((L * t).mid()) for t in ts]; rx = max(float((L * t).rad()) for t in ts)
        rows = sph_j_rows(xs, N, p)
        Bm = arb_mat(Np, n); BW = arb_mat(n, Np)
        wts = []
        for k, (r, w) in enumerate(nodes):
            t = ts[k]; wk = h * w / arb.pi(); gk = a_of_t(t) - g
            wts.append((wk, gk))
        for ii, d in enumerate(degs):
            sgn = -1 if ((d - p) // 2) % 2 else 1
            pre = sgn * (2 * L * (2 * d + 1)).sqrt()
            rr = rows[d]
            for k in range(n):
                val = pre * (rr[k] + arb(0, rx))      # |j_d'| <= 1 on R
                Bm[ii, k] = val; BW[k, ii] = val * wts[k][0] * wts[k][1]
        Ja += Bm * BW
        for ii in range(Np):
            bb[ii] += sum((Bm[ii, k] * Bm[0, k] * wts[k][0] for k in range(n)), arb(0))
        info.append(dict(panel=i, t0=float(ta), t1=float(tb), n=n, u=u, gmax=float(gmax.upper()), err_entry=float(eJ.upper())))
    save(path, dict(Ja=mat_to_fixed(Ja), bb=[(int(x.mid().man_exp()[0]), int(x.mid().man_exp()[1]), float(x.rad())) for x in bb],
                    errJ_entry=errJ, errb_entry=errb, info=info, seconds=time.time() - t0))
    return gid, 'done %.0fs' % (time.time() - t0)

# ---------------------------------------------------------------- assembly of T
def opnorm_bound(M):
    """||M||_2 <= ||M||_F, rigorous upper (from |entries| upper bounds)."""
    n, k = M.nrows(), M.ncols(); s = arb(0)
    for i in range(n):
        for j in range(k):
            s += arb(abs(M[i, j]).upper()) ** 2
    return arb(s.sqrt().upper())

def assemble(P, a, outdir, cst):
    ctx.prec = WP
    L, N, p, g, m = P['L'], P['N'], P['p'], P['g'], P['m']; degs = list(range(p, N, 2)); Np = len(degs)
    cells = sorted(f for f in os.listdir(os.path.join(outdir, 'cells')) if f.endswith('.pkl'))
    Cp = None; Cp2 = None
    for f in cells:
        D = pickle.load(open(os.path.join(outdir, 'cells', f), 'rb'))
        a1 = fixed_to_mat(D['Cp'], D['extra_Cp']); a2 = fixed_to_mat(D['Cp2'], D['extra_Cp2'])
        Cp = a1 if Cp is None else Cp + a1; Cp2 = a2 if Cp2 is None else Cp2 + a2
    Cp = (Cp + Cp.transpose()) / 2; Cp2 = (Cp2 + Cp2.transpose()) / 2
    Ja = None; bb = [arb(0)] * Np; errJ = 0.0; errb = 0.0
    for f in sorted(x for x in os.listdir(os.path.join(outdir, 'panels')) if x.endswith('.pkl')):
        D = pickle.load(open(os.path.join(outdir, 'panels', f), 'rb'))
        a1 = fixed_to_mat(D['Ja']); Ja = a1 if Ja is None else Ja + a1
        bb = [bb[i] + arb(man) * arb(2) ** ex + arb(0, r) for i, (man, ex, r) in enumerate(D['bb'])]
        errJ += D['errJ_entry']; errb += D['errb_entry']
    Ja = (Ja + Ja.transpose()) / 2
    I = arb_mat([[arb(1) if i == j else arb(0) for j in range(Np)] for i in range(Np)])
    A = g * I - Cp; B = g * g * I - 2 * g * Cp + Cp2
    # pole vector: int v_d e^{x/2} = sqrt(2L(2d+1)) i_d(L/2), i_d(z) = sqrt(pi/(2z)) I_{d+1/2}(z)
    z = L / 2
    pv = [(2 * L * (2 * d + 1)).sqrt() * (arb.pi() / (2 * z)).sqrt() * z.bessel_i(arb(d) + arb(1) / 2) for d in degs]
    psign = 2 if p == 0 else -2
    delta = arb(2) ** (-a.delta_log2)
    hv = [((1 - delta) if i == 0 else arb(0)) - bb[i] for i in range(Np)]
    uu = arb(a.u0) if p == 0 else arb(a.u1)
    J = Ja + arb_mat([[psign * pv[i] * pv[j] + uu * hv[i] * hv[j] for j in range(Np)] for i in range(Np)])
    # quadrature budget in operator norm: ||dJband|| <= Np * errJ ; U: ||u(hh'-h0h0')|| <= u(2||h|| eh + eh^2), eh = sqrt(Np) errb
    hn = math.sqrt(sum(float(abs(x).upper()) ** 2 for x in hv))
    eh = math.sqrt(Np) * errb
    eps_J = Np * errJ + float(uu) * (2 * hn * eh + eh * eh)
    X = (B - m * A).mid().solve((A - m * I).mid(), algorithm='approx').mid(); XT = X.transpose()
    W = X + XT - XT * A * X + (I - A * X - XT * A + XT * B * X) / m
    W = (W + W.transpose()) / 2
    Y = W.mid().solve(I, algorithm='approx').mid(); Y = ((Y + Y.transpose()) / 2).mid()
    R = I - Y * W
    nR = opnorm_bound(R); nY = opnorm_bound(Y); assert bool(nR < 1)
    eps_inv = float((nR * nY / (1 - nR)).upper())
    T = Y + J - cst['alpha'] * I - cst['beta'] * (B - A * A)
    T = (T + T.transpose()) / 2
    maxradT = max(float(T[i, j].rad()) for i in range(Np) for j in range(Np))
    # Guard (root cause of run 1's crash): a non-finite ball anywhere in J (Arb bessel_j at the turning point) makes every
    # entry of T = [0 +/- inf]; mid(T) = 0 and the later solve dies with "singular matrix".  Refuse to save such a T.
    if not (math.isfinite(maxradT) and math.isfinite(eps_J) and math.isfinite(eps_inv)):
        raise FloatingPointError('non-finite T: maxradT=%r eps_J=%r eps_inv=%r (check J panels / bessel top values)' % (maxradT, eps_J, eps_inv))
    budget = dict(eps_J=eps_J, eps_inv=eps_inv, normR=float(nR.upper()), normY=float(nY.upper()), maxradT=maxradT,
                  errJ_entry=errJ, errb_entry=errb, Np=Np)
    save(os.path.join(outdir, 'T.pkl'), dict(T=mat_to_fixed(T), budget=budget))
    return T, budget

# ---------------------------------------------------------------- certificate
def submat(M, r0, r1, c0, c1):
    return arb_mat([[M[i, j] for j in range(c0, c1)] for i in range(r0, r1)])

def ldl_congruence(M, prec):
    """Mid-precision recursive block LDL^T.  Returns (Z, d): Z exact unit lower triangular (arb_mat of exact midpoints),
    with Z M Z^T ~ diag(d).  Only an approximation; rigour comes from the ball check afterwards."""
    old = ctx.prec; ctx.prec = prec
    try:
        n = M.nrows()
        if n <= 48:
            A = [[M[i, j].mid() for j in range(n)] for i in range(n)]
            Lm = [[arb(0)] * n for _ in range(n)]; d = [arb(0)] * n
            for j in range(n):
                s = A[j][j] - sum((Lm[j][k] * Lm[j][k] * d[k] for k in range(j)), arb(0)); d[j] = s.mid()
                Lm[j][j] = arb(1)
                for i in range(j + 1, n):
                    t = A[i][j] - sum((Lm[i][k] * Lm[j][k] * d[k] for k in range(j)), arb(0))
                    Lm[i][j] = (t / d[j]).mid()
            Lmat = arb_mat(Lm)
            # Z = L^{-1} (unit lower), forward substitution in mid arithmetic
            Z = [[arb(0)] * n for _ in range(n)]
            for i in range(n):
                Z[i][i] = arb(1)
                for j in range(i):
                    Z[i][j] = (-sum((Lm[i][k] * Z[k][j] for k in range(j, i)), arb(0))).mid()
            return arb_mat(Z), d
        n1 = n // 2
        M11 = submat(M, 0, n1, 0, n1); M21 = submat(M, n1, n, 0, n1); M22 = submat(M, n1, n, n1, n)
        Z1, d1 = ldl_congruence(M11, prec)
        Dinv = arb_mat([[(1 / d1[i]).mid() if i == j else arb(0) for j in range(n1)] for i in range(n1)])
        Kf = (M21 * Z1.transpose() * Dinv * Z1).mid()          # M21 M11^{-1}
        S = (M22 - Kf * M21.transpose()).mid()
        Z2, d2 = ldl_congruence(S, prec)
        Lo = (-(Z2 * Kf)).mid()
        Z = [[arb(0)] * n for _ in range(n)]
        for i in range(n1):
            for j in range(n1): Z[i][j] = Z1[i, j]
        for i in range(n - n1):
            for j in range(n1): Z[n1 + i][j] = Lo[i, j]
            for j in range(n - n1): Z[n1 + i][n1 + j] = Z2[i, j]
        return arb_mat(Z), d1 + d2
    finally:
        ctx.prec = old

def certify(T, ell, prec, shift=None):
    """Returns dict: pd (bool, rigorous: diag dominance of S = Z(T - shift I)Z^T), neg_index (rigorous negative diagonal), etc.
    shift: exact arb point actually subtracted (defaults to arb(ell)); ell is kept only as a label."""
    old = ctx.prec; ctx.prec = prec
    try:
        n = T.nrows()
        I = arb_mat([[arb(1) if i == j else arb(0) for j in range(n)] for i in range(n)])
        sh = arb(ell) if shift is None else shift
        Ts = T - sh * I
        t0 = time.time(); Z, d = ldl_congruence(Ts, prec); t1 = time.time()
        S = Z * Ts * Z.transpose(); t2 = time.time()
        neg = [i for i in range(n) if bool(S[i, i] < 0)]
        dom_ok = True; worst = None; minpiv = None
        for i in range(n):
            off = sum((arb(abs(S[i, j]).upper()) for j in range(n) if j != i), arb(0))
            marg = S[i, i] - off
            if not bool(marg > 0): dom_ok = False
            r = float(S[i, i].mid()); minpiv = r if minpiv is None or r < minpiv else minpiv
            q = float((off / abs(S[i, i])).upper()) if bool(S[i, i] != 0) else float('inf')
            worst = q if worst is None or q > worst else worst
        return dict(ell=ell, shift=sh.str(20, radius=False), shift_minus_ell=(float((sh - arb(Fraction(str(ell)).numerator) / arb(Fraction(str(ell)).denominator)).mid()) if shift is not None else 0.0), pd_certified=dom_ok, neg_pivots=neg[:10], n_neg=len(neg), min_pivot=minpiv, worst_offdiag_ratio=worst,
                    t_ldl=t1 - t0, t_ball=t2 - t1)
    finally:
        ctx.prec = old

def lam_estimate(T, iters=30):
    """NUMERICAL: smallest eigenvalue estimate by inverse iteration on mid(T) (assumes the smallest is nearest 0)."""
    n = T.nrows(); Tm = T.mid()
    I = arb_mat([[arb(1) if i == j else arb(0) for j in range(n)] for i in range(n)])
    Ti = Tm.solve(I, algorithm='approx').mid()
    v = arb_mat([[arb(1) / (i + 1)] for i in range(n)]); lam = None
    for _ in range(iters):
        u = (Ti * v).mid()
        num = sum((u[i, 0] * v[i, 0] for i in range(n)), arb(0)); den = sum((u[i, 0] ** 2 for i in range(n)), arb(0))
        lam = num / den
        nu = den.sqrt(); v = arb_mat([[(u[i, 0] / nu).mid()] for i in range(n)])
    return float(lam.mid())

# ---------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['consts', 'cells', 'panels', 'assemble', 'certify', 'all'])
    ap.add_argument('--L', default='71/50'); ap.add_argument('--Omega', default='1000'); ap.add_argument('--s', default='763/200')
    ap.add_argument('--m', default='1/2'); ap.add_argument('--N', type=int, default=2376); ap.add_argument('--parity', type=int, default=0)
    ap.add_argument('--u0', type=int, default=607); ap.add_argument('--u1', type=int, default=200); ap.add_argument('--delta_log2', type=int, default=4200)
    ap.add_argument('--workers', type=int, default=3); ap.add_argument('--chunk', type=int, default=96)
    ap.add_argument('--out', required=True); ap.add_argument('--ells', default='')
    ap.add_argument('--cert_prec', type=int, default=576)
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out, 'cells'), exist_ok=True); os.makedirs(os.path.join(a.out, 'panels'), exist_ok=True)
    P = setup(a)
    cst, checks = consts(P, a)
    log = lambda *x: print(time.strftime('%H:%M:%S'), *x, flush=True)
    if a.stage in ('consts', 'all'):
        log('consts', json.dumps({k: v.str(8) for k, v in cst.items()})); log('checks', checks)
        json.dump(dict(consts={k: v.str(12) for k, v in cst.items()}, checks=checks, args=vars(a)), open(os.path.join(a.out, 'consts.json'), 'w'), indent=1)
    import multiprocessing as mpc
    if a.stage in ('cells', 'all'):
        bpa = breakpoints(P)
        log('cells', len(bpa) - 1, [round(float(x), 6) for x in bpa])
        GL.update(P=P, chunk=a.chunk, out=a.out, bpa=bpa); jobs = list(range(len(bpa) - 1))
        with mpc.get_context('fork').Pool(a.workers, maxtasksperchild=1) as pool:
            for r in pool.imap_unordered(cell_job, jobs): log('cell', r)
    if a.stage in ('panels', 'all'):
        plan = panel_plan(P)
        log('panels', len(plan), 'nodes', sum(x[3] for x in plan))
        G = max(1, a.workers * 2); groups = [plan[i::G] for i in range(G)]
        GL.update(P=P, out=a.out, groups=groups); jobs = list(range(G))
        with mpc.get_context('fork').Pool(a.workers, maxtasksperchild=1) as pool:
            for r in pool.imap_unordered(panel_job, jobs): log('panel group', r)
    if a.stage in ('assemble', 'all'):
        T, budget = assemble(P, a, a.out, cst); log('assembled', budget)
    if a.stage in ('certify', 'all'):
        D = pickle.load(open(os.path.join(a.out, 'T.pkl'), 'rb')); budget = D['budget']
        if not math.isfinite(D['T']['rad']) or not all(math.isfinite(budget[k]) for k in ('eps_J', 'eps_inv', 'maxradT')):
            raise FloatingPointError('T.pkl has non-finite radii/budget %r; re-run panels+assemble' % budget)
        T = fixed_to_mat(D['T'])
        sha = hashlib.sha256(json.dumps(D['T']['ints'][:1000] + [len(D['T']['ints'])]).encode()).hexdigest()
        res = dict(budget=budget, args=vars(a))
        est = lam_estimate(T); res['lam_est_numerical'] = est; log('lam_est (NUMERICAL, inverse iteration)', est)
        ells = a.ells.split(',') if a.ells else [repr(0.999 * est), repr(1.005 * est)]
        # Rigorous budget: T_true = W^{-1} + J_exact - alpha - beta(B - A^2) differs from the stored ball matrix T by at most
        # eps_inv (approximate inverse, Neumann bound) + eps_J (Gauss-Legendre remainder of the J band and of h, Bernstein
        # ellipses) in operator norm.  The floats were formed with upward-safe factors; doubling covers double rounding.
        extra = budget['eps_J'] + budget['eps_inv']
        res['extra_scalar_budget'] = extra
        ctx.prec = a.cert_prec
        ext = arb(2 * extra)
        res['certs'] = []
        for es in ells:
            fr = Fraction(es.strip()); ell_b = arb(fr.numerator) / arb(fr.denominator)       # ball containing the decimal ell
            # PD claim:  T - s_up I > 0 with s_up >= ell + ext  =>  lambda_min(T_true) > ell   (CERTIFIED lower bound)
            s_up = (ell_b + ext).upper()
            r = certify(T, es.strip(), a.cert_prec, shift=s_up); r['mode'] = 'pd_claim: lambda_min(T_true) > ell iff pd_certified'
            res['certs'].append(r); log('certify', r)
            # FAIL claim: certified negative diagonal of Z(T - s_lo I)Z^T, s_lo <= ell - ext  =>  lambda_min(T_true) < ell
            s_lo = (ell_b - ext).lower()
            r = certify(T, es.strip(), a.cert_prec, shift=s_lo); r['mode'] = 'fail_claim: lambda_min(T_true) < ell iff n_neg > 0'
            res['certs'].append(r); log('certify', r)
        json.dump(res, open(os.path.join(a.out, 'certify_%s.json' % time.strftime('%H%M%S')), 'w'), indent=1)

if __name__ == '__main__':
    main()
