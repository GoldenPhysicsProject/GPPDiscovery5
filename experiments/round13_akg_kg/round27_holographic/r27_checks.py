"""Round 27 checks for parts (1), (4) and the conservation-law candidates. No zeta zeros are proof inputs; the
off-line zeros of the controls (from Round 25 meta.json) are used only to locate their fake resonances.
Labels are given in ROUND27.md."""
import json, math, mpmath as mp, numpy as np
from r22_models import lam_table, KAP, TH
mp.mp.dps = 30
out = {}
CHI_M4 = (0, 1, 0, -1); CHI_5 = (0, 1, -1, -1, 1)
CHI_M20 = [CHI_M4[n % 4] * CHI_5[n % 5] for n in range(20)]
DHc = [0, 1, KAP, -KAP, -1]
LAM = {
 'ZETA': lambda s: mp.pi**(-s/2) * mp.gamma(s/2) * mp.zeta(s),
 'FTH':  lambda s: mp.pi**(-s) * mp.gamma((s+TH)/2) * mp.gamma((s-TH)/2) * mp.zeta(s+TH) * mp.zeta(s-TH),
 'EPS':  lambda s: (mp.sqrt(20)/(2*mp.pi))**s * mp.gamma(s) * (mp.zeta(s)*mp.dirichlet(s, CHI_M20) + mp.dirichlet(s, CHI_M4)*mp.dirichlet(s, CHI_5)),
 'DH':   lambda s: (5/mp.pi)**((s+1)/2) * mp.gamma((s+1)/2) * mp.dirichlet(s, DHc),
}
# (0) functional equation + reality of each completed function (the Klein symmetry, all four models)
for k, F in LAM.items():
    s0 = mp.mpc(0.31, 17.3)
    out[f'FE_{k}'] = float(abs(F(s0) - F(1 - s0)) / abs(F(s0)))
    out[f'reality_{k}'] = float(abs(F(mp.conj(s0)) - mp.conj(F(s0))) / abs(F(s0)))
# (1a) fake modular scattering phi_M(s) = Lam_M(2s-1)/Lam_M(2s) is unitary on Re s = 1/2 for every model
for k, F in LAM.items():
    out[f'unitarity_{k}'] = max(float(abs(abs(F(2*(0.5+1j*t)-1)/F(2*(0.5+1j*t))) - 1)) for t in (0.7, 3.1, 9.4, 21.0))
# (1b) resonances s0 = rho/2 and the cusp L^2 exponent: residue ~ y^{1-s0}, int^inf |y^{1-s0}|^2 y^-2 dy < inf  iff  Re s0 > 1/2
meta = json.load(open('../round25_viz/meta.json'))
res = {'ZETA': [0.5, meta['z_im'][0]], 'FTH': [0.75, meta['z_im'][0]], 'EPS': meta['eps_off'][0], 'DH': meta['dh_off'][0]}
out['resonance'] = {}
for k, (b, g) in res.items():
    s0 = mp.mpc(b, g) / 2
    val = abs(LAM[k](2*s0))
    out['resonance'][k] = dict(rho=[b, g], s0=[float(b/2), float(g/2)], Re_s0_minus_quarter=float(b/2 - 0.25),
                               cusp_integrand_exponent=float(-b), L2_normalizable=bool(b > 1), check_Lam_at_rho=float(val))
# (4a) finite-prime truncation of the modular scattering matrix: phi_P(s) = Gamma-ratio * prod_{p<=P} (1-p^{-2s})/(1-p^{1-2s})
def gr(s): return mp.pi**(-s/2) * mp.gamma(s/2)
def phiP(s, P):
    v = gr(2*s-1)/gr(2*s)
    for p in mp.primerange(2, P+1) if hasattr(mp, 'primerange') else [q for q in range(2, P+1) if all(q % d for d in range(2, int(q**0.5)+1))]:
        v *= (1 - mp.mpf(p)**(-2*s)) / (1 - mp.mpf(p)**(1-2*s))
    return v
out['truncated_phi_modulus_on_line'] = {str(P): [float(abs(phiP(0.5+1j*t, P))) for t in (3.1, 9.4, 21.0)] for P in (2, 13, 101, 1009)}
out['truncated_phi_pole_on_line_p2'] = dict(t=float(mp.pi/mp.log(2)), abs_phi_near=float(abs(phiP(0.5+1j*(mp.pi/mp.log(2)+1e-8), 13))))
# (4b) local Tate gammas gamma_p(s) = L_p(1-s)/L_p(s): unitary on the line for each p (all four 'local' types)
def gam_loc(s, p, alphas): return mp.fprod([(1 - a*mp.mpf(p)**(-s)) / (1 - a*mp.mpf(p)**(s-1)) for a in alphas])
out['local_unitarity_zeta'] = max(float(abs(abs(gam_loc(0.5+1j*t, p, [1])) - 1)) for p in (2, 3, 101) for t in (1.3, 7.7))
out['local_unitarity_fth'] = max(float(abs(abs(gam_loc(0.5+1j*t, p, [mp.mpf(p)**TH, mp.mpf(p)**-TH])) - 1)) for p in (2, 3, 101) for t in (1.3, 7.7))
# (4c) the per-prime U(1) phase does not converge as P -> inf on the line: arg prod_{p<=P} gamma_p(1/2+it)
N = 2_000_000; sv = np.ones(N+1, bool); sv[:2] = False
for p in range(2, int(N**0.5)+1):
    if sv[p]: sv[p*p::p] = False
P = np.nonzero(sv)[0].astype(float)
t = 10.0
ph = np.angle((1 - P**(-0.5-1j*t)) / (1 - P**(-0.5+1j*t)))   # each factor has modulus 1
cum = np.cumsum(ph)
out['phase_partial_sums_t10'] = {str(int(X)): float(cum[np.searchsorted(P, X)-1]) for X in (1e2, 1e3, 1e4, 1e5, 1e6, 2e6)}
# deterministic drift: sum_p p^{-1/2-it} ~ li(X^{1/2-it}) (prime number theorem); phase ~ -2 Im li(X^{1/2-it}) grows like sqrt(X)/(t log X)
out['phase_drift_prediction_2|li(X^(1/2+it))|'] = {str(int(X)): float(2*abs(mp.li(mp.mpf(X)**(0.5+1j*t)))) for X in (1e2, 1e4, 1e6, 2e6)}
# (C) conservation-law candidates, per model
lz = lam_table('ZETA', 5000); tabs = {k: lam_table(k, 5000) for k in ('FTH', 'EPS', 'DH')}
isppow = lz > 0
cand = {}
for k in ('ZETA', 'FTH', 'EPS', 'DH'):
    lm = lz if k == 'ZETA' else tabs[k]
    ratio = lm[isppow] / lz[isppow]
    cand[k] = dict(
        support_only_prime_powers=bool(np.all(np.abs(lm[~isppow][2:]) < 1e-9)),
        Lambda_nonneg=bool(np.all(lm[2:] > -1e-9)),
        n_negative_Lambda_le5000=int(np.sum(lm[2:] < -1e-9)),
        charge_ratio_min=float(ratio.min()), charge_ratio_max=float(ratio.max()))
cand['ZETA']['satake_moduli'] = [1.0]; cand['FTH']['satake_moduli'] = [2**TH, 2**-TH]
cand['EPS']['satake_moduli'] = 'undefined (no Euler product)'; cand['DH']['satake_moduli'] = 'undefined (no Euler product)'
out['conservation_candidates'] = cand
# winding of each local gamma over one period of t (argument principle per prime): same (0) for zeta and F_theta
def wind(alphas, p, n=4000):
    ts = np.linspace(0, 2*np.pi/np.log(p), n)
    z = np.array([complex(gam_loc(0.5+1j*tt, p, alphas)) for tt in ts])
    return float(np.round(np.sum(np.diff(np.unwrap(np.angle(z))))/(2*np.pi), 6))
out['local_winding'] = dict(zeta_p2=wind([1], 2), fth_p2=wind([mp.mpf(2)**TH, mp.mpf(2)**-TH], 2))
# Lee-Yang per prime: zeros of 1/L_p in z = p^{-s}: zeta z=1 (|z|=1); F_theta z = p^{-+theta} (|z| != 1)
out['lee_yang_zero_moduli_p2'] = dict(ZETA=[1.0], FTH=[float(2**-TH), float(2**TH)])
print(json.dumps(out, indent=1))
json.dump(out, open('r27_checks.json', 'w'), indent=1)
