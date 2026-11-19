"""Round 16 Part 3: the Stepanov port tried on Z, tested on zeta and every fake.  Integer / own-gamma data only.
(a) Hasse remainder envelope: from Round-15 one-shift blocks (l = 0.6), |b_s(g)|/d(g) per unit s-window and its growth
    exponent theta_W (b_s normalisation: trivial/Chebyshev = 1/2, RH <=> 0); theta_psi = 1/2 + theta_W.
(b) Chebyshev's add+mult identity  sum_n c(n) floor(x/n) = log floor(x)!   (Lambda * 1 = log: unique factorisation in Z)
    and the binomial step  sum_n c(n)(floor(2x/n) - 2 floor(x/n)) = log C(2x, x) <= 2x log 2 ; residual per model.
(c) Dimension-count (large sieve) step: Montgomery-Vaughan  sum_{x<p<=x+y} log p <= 2y log(x+y)/log y ; max ratio per model.
(d) Landau step (Z-analogue of the twist / lower-bound step): one-sided bound + no real zero in (1/2,1) => two-sided.
    Real-axis values of each model's L-function on (1/2, 1) (mpmath, from its own Dirichlet series / characters).
(e) Tensor-power analogue: progressions s0, 2 s0, 3 s0 ... ; does |b_{s0}| <= d propagate to |b_{n s0}| <= d?  Hankel
    numerical rank (genus) and matrix-pencil growth rate e^{delta} of n -> b_{s0 + n h} (an OUTPUT, zeros never input)."""
import sys, os, json, math, numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); R15 = os.path.join(HERE, '..', 'round15_hasse')
sys.path.insert(0, R15); sys.path.insert(0, os.path.join(HERE, '..', 'round14b_euler')); sys.path.insert(0, os.path.join(HERE, '..'))
import hasse_core as hc
NAMES = ['ZETA', 'FTAU', 'EPS', 'DH', 'BEUR2', 'FTH']
out = {}

# (a) envelope
env = {}
for nm in NAMES:
    D = json.load(open(os.path.join(R15, f'blocks_{nm}_l0.6.json'))); rows = D['rows']
    rec = {}
    for tag in ('sine1', 'ground'):
        s = np.array([r['s'] for r in rows]); b = np.array([r['fixed'][tag]['b'] for r in rows]); d = rows[0]['fixed'][tag]['d']
        pole = np.array([r['fixed'][tag]['pole'] for r in rows])
        wins = []
        for k in range(1, 8):
            m = (s >= k) & (s < k + 1)
            if m.sum(): wins.append(dict(win=[k, k+1], max_b_over_d=float(np.abs(b[m]).max()/d), max_b_over_pole=float((np.abs(b[m])/np.where(pole[m] != 0, pole[m], np.nan)).max()) if np.any(pole[m]) else None))
        c = np.array([w['win'][0] + 0.5 for w in wins]); y = np.log([w['max_b_over_d'] for w in wins]); sel = c >= 2
        th = float(np.polyfit(c[sel], y[sel], 1)[0]) if sel.sum() >= 3 else None
        rec[tag] = dict(d=d, windows=wins, theta_W_fit=th, theta_psi_fit=(0.5 + th) if th is not None else None,
                        max_b_over_d=float(np.abs(b).max()/d))
    sr = np.array([r['s'] for r in rows]); rho = np.array([r['rho'] if r['rho'] is not None else np.nan for r in rows])
    rw = []
    for k in range(1, 8):
        m = (sr >= k) & (sr < k + 1)
        if m.sum(): rw.append((k + 0.5, float(np.nanmax(rho[m]))))
    cc = np.array([a for a, _ in rw]); yy = np.log([b_ for _, b_ in rw]); sel = cc >= 2
    rec['rho'] = dict(windows=rw, theta_W_fit=float(np.polyfit(cc[sel], yy[sel], 1)[0]), max_rho=float(np.nanmax(rho)))
    env[nm] = rec
    print('(a)', nm, 'rho-env slope=%.3f max rho=%.3f' % (rec['rho']['theta_W_fit'], rec['rho']['max_rho']), 'theta_W(sine1)=%.3f theta_W(ground)=%.3f max|b|/d=%.3f' % (rec['sine1']['theta_W_fit'], rec['ground']['theta_W_fit'], rec['sine1']['max_b_over_d']), flush=True)
out['a_envelope'] = env

# (b) Chebyshev identity + binomial step, (c) MV window step
models = {nm: hc.MODELS[nm]() for nm in NAMES}
cheb, mv = {}, {}
for nm, mo in models.items():
    n = mo['n']; c = mo['c']; res = []; binom = []
    for x in (100, 300, 1000, 2000, 2700):
        S = float(np.sum(c[n <= x]*np.floor(x/n[n <= x]))); L = math.lgamma(x + 1)
        res.append(dict(x=x, identity_residual_over_x=(S - L)/x))
        Bm = float(np.sum(c[n <= 2*x]*(np.floor(2*x/n[n <= 2*x]) - 2*np.floor(x/n[n <= 2*x]))))
        lb = math.lgamma(2*x + 1) - 2*math.lgamma(x + 1)
        psi = lambda X: float(np.sum(c[n <= X]))
        binom.append(dict(x=x, binom_residual_over_x=(Bm - lb)/x, psi2x_minus_psix_over_2xlog2=(psi(2*x) - psi(x))/(2*x*math.log(2))))
    cheb[nm] = dict(identity=res, binomial=binom)
    # MV: windows (x, x+y], y = x^0.5, x^0.75, and the Hasse-window width x(e^l - e^-l), l = 0.6
    pr = {}
    for lab, yf in (('y=x^0.5', lambda x: x**0.5), ('y=x^0.75', lambda x: x**0.75), ('y=x(e^.6-e^-.6)', lambda x: x*(math.exp(0.6) - math.exp(-0.6)))):
        worst = 0.0; wx = None
        for x in np.geomspace(50, n.max()/3.5, 300):
            y = yf(x); m = (n > x) & (n <= x + y); S = float(np.sum(c[m]))
            r = S/(2*y*math.log(x + y)/math.log(y)); 
            if r > worst: worst, wx = r, float(x)
        pr[lab] = dict(max_ratio=worst, at_x=wx)
    mv[nm] = pr
    print('(b,c)', nm, 'cheb id resid/x @2700 = %.4f' % res[-1]['identity_residual_over_x'], 'binom resid/x = %.4f' % binom[-1]['binom_residual_over_x'],
          '(psi2x-psix)/(2xlog2) = %.3f' % binom[-1]['psi2x_minus_psix_over_2xlog2'], 'MV max ratios', {k: round(v['max_ratio'], 3) for k, v in pr.items()}, flush=True)
out['b_chebyshev'] = cheb; out['c_MV'] = mv

# (d) real-axis values on (1/2, 1)
mp.mp.dps = 20
chi = lambda tab: tab
X4 = [0, 1, 0, -1]; X5 = [0, 1, -1, -1, 1]
X20 = [0] + [0]*19
for k in range(20): X20[k] = 0 if math.gcd(k, 20) > 1 else (X4[k % 4]*X5[k % 5])   # chi_{-20} = chi_{-4} chi_5
KAP = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
Xc = [0, 1, 1j, -1j, -1]
def val(nm, s):
    s = mp.mpf(s)
    if nm in ('ZETA',): return mp.zeta(s)
    if nm == 'EPS': return mp.zeta(s)*mp.dirichlet(s, X20) + mp.dirichlet(s, X4)*mp.dirichlet(s, X5)   # Z_{x^2+5y^2} (sum over (x,y)!=0)
    if nm == 'DH':
        L = mp.dirichlet(s, Xc); return mp.re(L) + KAP*mp.im(L)
    if nm == 'FTH': return mp.zeta(s - 0.2)*mp.zeta(s + 0.2)
    if nm == 'FTAU': return abs(mp.zeta(s + 3j))**2
    return None
# sanity: EPS formula vs direct lattice sum at s = 2
import harness as hh
rr = hh.epstein_model()['r']; direct = float(np.sum(rr[1:]*np.arange(1, len(rr), dtype=float)**-2.0))
real = {'EPS_formula_check_s2': dict(formula=float(val('EPS', 2)), lattice_sum_truncated=direct)}
for nm in ('ZETA', 'EPS', 'DH', 'FTH', 'FTAU'):
    grid = np.linspace(0.505, 0.995, 99); v = [float(val(nm, x)) for x in grid]
    sc = [float((grid[i] + grid[i+1])/2) for i in range(len(v) - 1) if v[i]*v[i+1] < 0]
    real[nm] = dict(min=min(v), max=max(v), sign_changes_near=sc, real_zero_in_half_one=bool(sc))
    print('(d)', nm, 'range [%.4g, %.4g]' % (min(v), max(v)), 'real zeros near', sc, flush=True)
real['BEUR2'] = 'N.A. (no continuation from finite data)'
out['d_real_axis'] = real

# (e) progressions and Hankel / matrix pencil
prog = {}
for nm in NAMES:
    D = json.load(open(os.path.join(R15, f'blocks_{nm}_l0.6.json'))); rows = [r for r in D['rows'] if not r['prime_shift']]
    s = np.array([r['s'] for r in rows]); b = np.array([r['fixed']['sine1']['b'] for r in rows]); d = rows[0]['fixed']['sine1']['d']
    u = b/d; rec = {}
    for s0 in (0.6, 0.8, 1.0, 1.2, 1.4):
        seq = []
        for k in range(1, 40):
            j = np.nonzero(np.abs(s - k*s0) < 1e-6)[0]
            if not len(j): break
            seq.append(float(u[j[0]]))
        first = next((k+1 for k, v in enumerate(seq) if abs(v) > 1), None)
        rec[f's0={s0}'] = dict(u=seq, holds_at_n1=abs(seq[0]) <= 1, first_n_violating=first)
    # Hankel numerical rank and matrix pencil on the uniform grid
    h = 0.02; K = len(u)//2; Hk = np.array([[u[i + j] for j in range(K)] for i in range(len(u) - K)])
    sv = np.linalg.svd(Hk, compute_uv=False); rk = int(np.sum(sv > 1e-6*sv[0]))
    r_use = min(rk, 40); Y0 = Hk[:, :-1]; Y1 = Hk[:, 1:]
    U, S_, Vt = np.linalg.svd(Y0, full_matrices=False); U, S_, Vt = U[:, :r_use], S_[:r_use], Vt[:r_use]
    z = np.linalg.eigvals(np.diag(1/S_) @ U.T @ Y1 @ Vt.T)
    # amplitude of each mode (least squares), keep modes carrying >= 1e-3 of the signal
    V = np.vander(z, len(u), increasing=True).T; amp = np.linalg.lstsq(V, u.astype(complex), rcond=None)[0]
    keep = np.abs(amp)*np.maximum(1, np.abs(z)**(len(u) - 1)) > 1e-3*np.abs(u).max()
    rate = np.log(np.abs(z[keep]))/h
    rec['hankel'] = dict(K=K, numerical_rank_rel1e6=rk, sv_head=[float(x) for x in sv[:6]/sv[0]])
    rec['pencil'] = dict(modes=int(keep.sum()), max_growth_delta=float(rate.max()), min_growth_delta=float(rate.min()),
                         n_modes_delta_gt_003=int(np.sum(rate > 0.03)), top_modes=[dict(delta=float(np.log(abs(zz))/h), freq=float(abs(np.angle(zz))/h), amp=float(abs(a)))
                                                                                   for zz, a in sorted(zip(z[keep], amp[keep]), key=lambda t: -np.log(abs(t[0])))[:6]])
    rec['pencil']['caveat'] = 'ill-conditioned: FTAU (no off-line zeros) shows a spurious delta; NOT used as evidence'
    rr = [r for r in D['rows'] if not r['prime_shift']]; sr = np.array([r['s'] for r in rr]); rh = np.array([r['rho'] for r in rr])
    for s0 in (0.6, 0.8, 1.0, 1.2, 1.4):
        seq = []
        for k in range(1, 40):
            j = np.nonzero(np.abs(sr - k*s0) < 1e-6)[0]
            if not len(j): break
            seq.append(float(rh[j[0]]))
        rec[f'rho_s0={s0}'] = dict(rho=seq, holds_at_n1=seq[0] <= 1, first_n_violating=next((k+1 for k, v in enumerate(seq) if v > 1), None))
    prog[nm] = rec
    print('(e)', nm, {k: (v['holds_at_n1'], v['first_n_violating']) for k, v in rec.items() if k.startswith('s0') or k.startswith('rho_')}, 'hankel rank', rk,
          'pencil max delta %.3f' % rec['pencil']['max_growth_delta'], flush=True)
# F_q comparison: Hankel rank of n -> a_n is exactly 2g
ff = json.load(open(os.path.join(HERE, 'stepanov_ff.json')))
hk = []
for cv in ff['curves']:
    om = []
    p, g = cv['p'], cv['g']
    # rebuild alphas from abs_omega is lossy; recompute a_n from stored N via formula is cheaper: use ext deviations sign-free -> use rank of synthetic
    hk.append(dict(curve=cv['curve']['tag'], genus=g, hankel_rank_expected=2*g))
out['e_progressions'] = prog
# (f) the PROVED Chebyshev-strength bound for zeta (binomial step only, no zeros):  |prime_s(g)| <= 2 log2 e^{3l/2} e^{s/2} ||g||^2
D = json.load(open(os.path.join(R15, 'blocks_ZETA_l0.6.json'))); l = 0.6; fb = []
for r in D['rows']:
    if r['prime_shift'] or abs(r['s']*2 - round(r['s']*2)) > 1e-9: continue
    fx = r['fixed']['sine1']; proven = abs(fx['pole']) + abs(fx['arch']) + 2*math.log(2)*math.exp(1.5*l)*math.exp(r['s']/2)
    fb.append(dict(s=r['s'], measured_abs_b=abs(fx['b']), d=fx['d'], proven_bound=proven, proven_over_d=proven/fx['d']))
out['f_proven_bound_zeta'] = fb
print('(f)', [(x['s'], round(x['proven_over_d'], 1), round(x['measured_abs_b']/x['d'], 3)) for x in fb[::2]]); out['e_ff_hankel_rank'] = hk
json.dump(out, open('zport.json', 'w'), indent=1, default=float)
print('DONE')
