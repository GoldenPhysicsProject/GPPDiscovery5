"""Round 33 analysis: the torus SOS on the Fejer (Kronecker-line) measure nu_L, for every model with its own mu.
f_c = sum_{n<=X} a_n chi_n - c chi_1,  a_n = mu_M(n) n^{-1/2},  c = 1/b_M(1/2)   (Arb constants from R32)
  <f_c,1>_{nu_L}   = m(L) - cL                      (the Moebius floor, R32; = m~ + intercept)
  ||f_c||^2_{nu_L} = G_c = 2 Q_c + L Mc(X)^2 - 2 Mc(X) m_c,   Q_c = int_0^L (M(e^u)-c)^2 du   (exact step integrals)
  Haar part H_c   = L * sum |coef|^2   (counted once);   remainder R_c = G_c - H_c
  SOS / Cauchy-Schwarz: (m - cL)^2 <= L Q_c <= L G_c.   psi side (ZETA): coefficients (Lambda-1) n^{-1/2}, c_Lam."""
import json, math, numpy as np, mpmath
K = json.load(open('../round32_mollified_floor/r32_constants.json'))
BZ = json.load(open('../round32_mollified_floor/r32_beur_zeta.json'))
mpmath.mp.dps = 30
c_lam = float(-mpmath.zeta(0.5, derivative=1)/mpmath.zeta(0.5) - mpmath.zeta(0.5))   # -zeta'/zeta(1/2) - zeta(1/2)
P = json.load(open('r33_pass_1e09.json')); C = json.load(open('r33_controls_1e07.json')); B = json.load(open('r33_beurling.json'))
def from_thresholds(D, key):
    T = np.array(D['T']); H = D['H']; Lg = np.arange(H, math.log(D['N']) + 1e-12, H)
    j = np.searchsorted(T, np.floor(np.exp(Lg)*(1 + 1e-13)).astype(np.int64)); NX = T[j].astype(float)
    g = lambda s: np.array(D['cum'][f'{key}_{s}'])[j]
    M = g('M'); l1 = np.log1p(1/NX); tail = Lg - np.log(NX)
    return Lg, g('dg'), g('q') - M*M*l1 + M*M*tail, g('r') - M*l1 + M*tail, M
models = {'ZETA': ('mu', P, K['ZETA']), 'ZETA2': ('f0.0', P, K['FTH0.0']), 'FTH0.1': ('f0.1', P, K['FTH0.1']),
          'FTH0.25': ('f0.25', P, K['FTH0.25']), 'EPS': ('EPS', C, K['EPS']), 'DH': ('DH', C, K['DH'])}
res, ser = {}, {}
def sos(name, L, dg, q, r, M, c, icpt):
    Qc = q - 2*c*r + c*c*L; mc = r - c*L; Mc = M - c
    Gc = 2*Qc + L*Mc**2 - 2*Mc*mc; Hc = L*(dg - 1 + (1 - c)**2); Rc = Gc - Hc
    eta = mc**2/(L*Qc); mt = mc - icpt
    s = dict(L=L, Qc=Qc, Gc=Gc, Hc=Hc, Rc=Rc, eta=eta, mt=mt, mc=mc); ser[name] = s
    sel = L >= 3
    def at(x):
        i = int(np.argmin(abs(L - x))); return dict(L=float(L[i]), Qc=float(Qc[i]), Gc=float(Gc[i]), Hc=float(Hc[i]), Rc=float(Rc[i]),
                                                  Rc_over_Hc=float(Rc[i]/Hc[i]), Qc_over_L=float(Qc[i]/L[i]), eta=float(eta[i]))
    pts = [at(x) for x in (6, 9, 12, 15, 18, float(L[-1]))]
    # local slope of Qc on the last third of the range and an exponential-rate fit of log Qc
    hi = L >= L[-1]*2/3
    slope = float(np.polyfit(L[hi], Qc[hi], 1)[0]); rate = float(np.polyfit(L[hi], np.log(np.maximum(Qc[hi], 1e-300)), 1)[0])
    res[name] = dict(c=c, intercept=icpt, Lmax=float(L[-1]), points=pts, Qc_slope_last_third=slope, logQc_rate_last_third=rate,
                     Rc_sign_L_ge_3=dict(min=float(Rc[sel].min()), max=float(Rc[sel].max()), frac_negative=float((Rc[sel] < 0).mean())),
                     CS_max_eta=float(eta[sel].max()), CS_holds=bool((eta[sel] <= 1 + 1e-12).all()) and bool((L*Gc >= mc**2 - 1e-9)[sel].all()),
                     mt_range=[float(mt[sel].min()), float(mt[sel].max())])
for name, (key, D, k) in models.items():
    L, dg, q, r, M = from_thresholds(D, key); sos(name, L, dg, q, r, M, k['slope'], k['intercept'])
L = np.array(B['L']); sos('BEUR', L, np.array(B['dg']), np.array(B['q']), np.array(B['r']), np.array(B['M']), BZ['slope'], BZ['intercept'])
# psi side for ZETA: vector (Lambda-1) n^{-1/2} - c_Lam chi_1
L, dg, q, r, M = from_thresholds(P, 'lam'); sos('ZETA_PSI', L, dg, q, r, M, c_lam, 0.0)
res['ZETA_PSI']['c_lam'] = c_lam
# tie to the R31 floor D (R32 series): D + <f,1> should be a constant kappa0 (+ O(e^{-L/2}))
z = np.load('../round32_mollified_floor/r32_series.npz'); Lz = z['ZETA_Lg']; Dz = z['ZETA_D']
s = ser['ZETA_PSI']; cr = np.interp(Lz, s['L'], s['mc']); k0 = Dz + cr; sel = Lz >= 8
res['ZETA_PSI']['D_plus_cross'] = dict(mean=float(k0[sel].mean()), spread=float(k0[sel].max() - k0[sel].min()))
# consistency with R32 m~
for name, nk in (('ZETA', 'ZETA'), ('ZETA2', 'ZETA2'), ('FTH0.1', 'FTH0.1'), ('FTH0.25', 'FTH0.25')):
    s = ser[name]; Lr = z[f'{nk}_Lg']; mr = z[f'{nk}_m1t']; res[name]['max_diff_vs_R32_mt'] = float(np.abs(np.interp(Lr, s['L'], s['mt']) - mr)[Lr >= 3].max())
json.dump(res, open('r33_results.json', 'w'), indent=1)
np.savez_compressed('r33_series.npz', **{f'{n}_{k}': v for n, s in ser.items() for k, v in s.items()})
for n, v in res.items():
    print(n, 'c=%.4f' % v['c'], 'Qc/L slope %.4g rate %.4g' % (v['Qc_slope_last_third'], v['logQc_rate_last_third']),
          'R<0 frac %.3f' % v['Rc_sign_L_ge_3']['frac_negative'], 'etamax %.3g' % v['CS_max_eta'], v['CS_holds'], v.get('max_diff_vs_R32_mt'))
    for p in v['points']: print('   ', {k: round(x, 5) for k, x in p.items()})
print(res['ZETA_PSI']['D_plus_cross'])
