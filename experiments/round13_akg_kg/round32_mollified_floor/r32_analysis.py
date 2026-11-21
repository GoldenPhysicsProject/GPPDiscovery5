"""Round 32 analysis: the R31 floor D_M and the Moebius floor m1~_M for ZETA, F_theta (0.1, 0.25), EPS, DH, BEUR.
  m1(L)  = sum_{d<=e^L} mu_M(d) d^{-1/2}(L - log d)          (Laplace transform 1/(b_M(s+1/2) s^2))
  m1~(L) = m1(L) - L/b(1/2) + b'(1/2)/b(1/2)^2               (Arb-certified constants, r32_constants.json)
  BEUR: no FE / closed form for b(1/2): its (slope, intercept) are least-squares FITTED on L in [6, Lmax] (labelled).
  D_M(L) = P_M(L) - S_M(L) (R31 floor; P from r32_common).  ZETA decomposition (identity A3 with A(Y) asymptotics):
  D = 4 sqrt(X)(1 - m^(1)) - 2 zeta(1/2) m2~ - zeta'(1/2) m1~ - E_mu + R0(L)."""
import json, math, sys, numpy as np
from r32_common import P_zeta, P_model, P_beur, c1, C2

PASS = sys.argv[1] if len(sys.argv) > 1 else 'r32_pass_1e09.json'
pz = json.load(open(PASS)); pc = json.load(open('r32_controls_1e07.json'))
pb = json.load(open('r32_beurling.json')); K = json.load(open('r32_constants.json'))
pz0 = json.load(open('r32_pass_1e09_th0.0.json')); KB = json.load(open('r32_beur_zeta.json'))
H = 0.005
out = {}

def series(src, pre, preP, Nmax, const, th=None, name=None):
    T = np.array(src['T']); cum = src['cum']
    A0, A1, A2 = (np.array(cum[f'{pre}_A{j}']) for j in range(3))
    P0, P1 = np.array(cum[f'{preP}_P0']), np.array(cum[f'{preP}_P1'])
    # L-grid
    Lg = np.arange(H, math.log(Nmax) + 1e-12, H)
    ig = np.searchsorted(T, np.floor(np.exp(Lg)*(1 + 1e-13)).astype(np.int64))
    ok = (ig < len(T)) & (T[np.minimum(ig, len(T)-1)] == np.floor(np.exp(Lg)*(1 + 1e-13)))
    Lg, ig = Lg[ok], ig[ok]
    m1 = Lg*A0[ig] - A1[ig]; m2 = 0.5*(Lg**2*A0[ig] - 2*Lg*A1[ig] + A2[ig])
    m1t = m1 - const['slope']*Lg - const['intercept']
    S = Lg*P0[ig] - P1[ig]
    Pp = np.array([P_zeta(L) if name == 'ZETA' else P_model(name, L, th) for L in Lg])
    D = Pp - S
    # integer samples k^4
    kmax = src['kmax']; ks = np.arange(2, kmax + 1); ik = np.searchsorted(T, ks.astype(np.int64)**4)
    Lk = 4*np.log(ks)
    m1k = Lk*A0[ik] - A1[ik] - const['slope']*Lk - const['intercept']
    Dk = np.array([P_zeta(L) if name == 'ZETA' else P_model(name, L, th) for L in Lk]) - (Lk*P0[ik] - P1[ik])
    return dict(Lg=Lg, m1t=m1t, m1=m1, m2=m2, D=D, ig=ig, ks=ks, Lk=Lk, m1k=m1k, Dk=Dk)

def windows(ks, v):
    res = []
    for Kk in range(2, ks[-1] + 1):
        K2 = math.ceil(Kk**1.5)
        if K2 > ks[-1]: break
        sel = (ks >= Kk) & (ks <= K2); d = v[sel]
        res.append(dict(K=Kk, K32=K2, min=float(d.min()), max=float(d.max()), maxabs=float(np.abs(d).max())))
    return res

def growth(Lg, v, L0=6.0):
    """window envelope W(L) = max_{[L, 1.5L]} |v| on the grid; slope of log W vs L for L >= L0."""
    Ls, W = [], []
    for L in np.arange(L0, Lg[-1]/1.5, 0.25):
        sel = (Lg >= L) & (Lg <= 1.5*L); Ls.append(L); W.append(float(np.abs(v[sel]).max()))
    Ls, W = np.array(Ls), np.array(W)
    sl = float(np.polyfit(Ls, np.log(W), 1)[0]) if len(Ls) > 3 else float('nan')
    return dict(L=Ls.tolist(), W=W.tolist(), slope=sl)

specs = [('ZETA', pz, 'mu', 'lam', pz['N'], K['ZETA'], None),
         ('ZETA2', pz0, 'f0.0', 'lf0.0', pz0['N'], K['FTH0.0'], 0.0),
         ('FTH0.1', pz, 'f0.1', 'lf0.1', pz['N'], K['FTH0.1'], 0.1),
         ('FTH0.25', pz, 'f0.25', 'lf0.25', pz['N'], K['FTH0.25'], 0.25),
         ('EPS', pc, 'EPS', 'EPS', pc['N'], K['EPS'], None),
         ('DH', pc, 'DH', 'DH', pc['N'], K['DH'], None)]
R = {}
for key, src, pre, preP, Nm, const, th in specs:
    name = 'FTH' if (key.startswith('FTH') or key == 'ZETA2') else key
    r = series(src, pre, preP, Nm, const, th, name); R[key] = r
    sel = r['Lg'] >= 3
    out[key] = dict(N=Nm, D_min_int=float(r['Dk'].min()), D_min_int_k=int(r['ks'][np.argmin(r['Dk'])]),
                    D_min_grid=float(r['D'][sel].min()), D_min_grid_L=float(r['Lg'][sel][np.argmin(r['D'][sel])]),
                    m1t_min_grid=float(r['m1t'][sel].min()), m1t_max_grid=float(r['m1t'][sel].max()),
                    m1t_min_int=float(r['m1k'][r['Lk'] >= 3].min()), m1t_max_int=float(r['m1k'][r['Lk'] >= 3].max()),
                    m1t_end=float(r['m1t'][-1]),
                    windows_m1t=windows(r['ks'], r['m1k']), windows_D=windows(r['ks'], r['Dk']),
                    growth_m1t=growth(r['Lg'], r['m1t']), growth_D=growth(r['Lg'], np.minimum(r['D'], 0) - 1e-300),
                    dyadic_max_mu=(src['dyadic_max'].get({'ZETA': 'mu', 'ZETA2': 'f0.0', 'FTH0.1': 'f0.1', 'FTH0.25': 'f0.25'}.get(key, key)))[:int(math.log2(Nm)) + 1])
    print(key, {k: out[key][k] for k in ('D_min_int', 'D_min_int_k', 'm1t_min_grid', 'm1t_max_grid', 'm1t_min_int', 'm1t_max_int')},
          'growth slope', round(out[key]['growth_m1t']['slope'], 3), flush=True)

# BEUR
Lb = np.array(pb['L']); A0, A1 = np.array(pb['A0']), np.array(pb['A1']); P0, P1 = np.array(pb['P0']), np.array(pb['P1'])
selb = (Lb >= H) & (Lb <= math.log(pb['NMAX']))
Lb, A0, A1, P0, P1 = Lb[selb], A0[selb], A1[selb], P0[selb], P1[selb]
m1b = Lb*A0 - A1
fit = Lb >= 6
cb = np.polyfit(Lb[fit], m1b[fit], 1)                       # reference only
m1bt = m1b - KB['slope']*Lb - KB['intercept']              # constants from the continuation (r32_beur_zeta.py)
cD = np.polyfit(Lb[fit], (np.array([P_beur(L) for L in Lb]) - (Lb*P0 - P1))[fit], 1)
Db = np.array([P_beur(L) for L in Lb]) - (Lb*P0 - P1) - np.polyval(cD, Lb)   # linear part FITTED (no FE)
R['BEUR'] = dict(Lg=Lb, m1t=m1bt, D=Db)
s3 = Lb >= 3.5
out['BEUR'] = dict(N=pb['NMAX'], fitted_slope=float(cb[0]), fitted_intercept=float(cb[1]), slope=KB['slope'], intercept=KB['intercept'],
                   D_linear_fit=cD.tolist(), m1t_at={str(l): float(m1bt[np.argmin(abs(Lb - l))]) for l in (6, 8, 10, 12, 14, 16, 18, 18.4)},
                   m1t_min_grid=float(m1bt[s3].min()), m1t_max_grid=float(m1bt[s3].max()),
                   D_min_grid=float(Db[s3].min()), D_max_grid=float(Db[s3].max()), growth_m1t=growth(Lb, m1bt),
                   nprimes=pb['nprimes'], PiB_minus_PiC=pb['PiB_minus_PiC'], dfs=pb['dfs'])
print('BEUR', {k: out['BEUR'][k] for k in ('fitted_slope', 'm1t_min_grid', 'm1t_max_grid', 'D_min_grid')}, flush=True)

# ZETA decomposition
r = R['ZETA']; z = K['ZETA']; z0, z1, inv2 = z['z0'], z['z1'], z['inv2']
T = np.array(pz['T']); B0, B1 = np.array(pz['cum']['mu_B0']), np.array(pz['cum']['mu_B1'])
L = r['Lg']; X = np.exp(L)
mone = L*B0[r['ig']] - B1[r['ig']]
T1 = 4*np.sqrt(X)*(1 - mone)
m2t = r['m2'] - 0.5*(L**2/z0 - 2*L*z1/z0**2 + inv2)
T2 = -2*z0*m2t; T3 = -z1*r['m1t']
_b = 2*np.arange(0, 400) + 0.5
R0 = np.array([-8 + C2 - z0*inv2 + z1**2/z0**2 + 4*math.exp(-l/2) - float(np.sum(np.exp(-_b*l)/_b**2)) for l in L])
Emu = T1 + T2 + T3 + R0 - r['D']
R['ZETA'].update(T1=T1, T2=T2, T3=T3, Emu=Emu)
s = L >= 4
out['ZETA_decomposition'] = {nm: dict(min=float(v[s].min()), max=float(v[s].max())) for nm, v in
                             (('4sqrtX(1-m^(1))', T1), ('-2z0 m2~', T2), ('-z1 m1~', T3), ('E_mu (residual)', Emu), ('D', r['D']))}
out['R0_const'] = float(-8 + C2 - z0*inv2 + z1**2/z0**2)
print('decomp', out['ZETA_decomposition'], flush=True)
json.dump(out, open('r32_results.json', 'w'))
np.savez_compressed('r32_series.npz', **{f'{k}_{f}': np.asarray(v) for k, d in R.items() for f, v in d.items()})
