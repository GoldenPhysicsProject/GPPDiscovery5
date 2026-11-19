"""Round 17B: the massless spin-2 (graviton) conformal primary version of "s IS Delta".

Map (PROVED by change of variables, CHECKED here): for each model M with Round-17 law Phi_M,
    graviton radial (energy) profile   a_M(omega) = 2 f_M(omega^2) = 2 omega^{-1} Phi_M(2 log omega)
    PSS energy-Mellin                  A_M(Delta) = int_0^oo omega^{Delta-1} a_M(omega) d omega = xi_M(Delta/2).
KG norm of the spin-2 wavepacket (per unit angular profile) ~ int |a|^2 omega d omega = 2 int Phi^2,
and Plancherel on PSS's principal series Re Delta = 1:  int |a|^2 omega d omega = (1/2pi) int |A(1+i nu)|^2 d nu.

Integer data only: a_M is built from theta / Dirichlet coefficients (laws17.Model). Closed forms of xi_M
(mpmath zeta / Dirichlet series / gamma) are an independent check. Zeros are NOT inputs; the one EPS zero
below is located by Newton from Round 17's printed label only to illustrate the GAP (sec. F).
"""
import json, math, time, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
mp.mp.dps = 25
from laws17 import Model, U_MAX, kappa

OUT = {}
t0 = time.time()

def chi_list(mod, f):
    return [mp.mpf(f(n)) for n in range(mod)]

CHI_M4 = chi_list(4, lambda n: [0, 1, 0, -1][n % 4])
CHI_5 = chi_list(5, lambda n: [0, 1, -1, -1, 1][n % 5])
CHI_M20 = chi_list(20, lambda n: CHI_M4[n % 4]*CHI_5[n % 5])

def xi_closed(name, s):
    s = mp.mpc(s)
    if min(abs(s), abs(s-1)) < mp.mpf('1e-12'):      # removable points: symmetric average
        h = mp.mpf('1e-10')
        return (xi_closed(name, s+h) + xi_closed(name, s-h))/2
    if name == 'ZETA':
        return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
    if name == 'EPS':   # Z_Q = zeta_K + L(chi_-4) L(chi_5), K = Q(sqrt -5), Q = m^2 + 5 n^2
        ZQ = mp.zeta(s)*mp.dirichlet(s, CHI_M20) + mp.dirichlet(s, CHI_M4)*mp.dirichlet(s, CHI_5)
        return s*(s-1)*(mp.sqrt(20)/(2*mp.pi))**s*mp.gamma(s)*ZQ
    if name == 'DH':
        k = kappa()
        return (5/mp.pi)**((s+1)/2)*mp.gamma((s+1)/2)*mp.dirichlet(s, [0, 1, k, -k, -1])
    raise ValueError(name)

def theta_tail(x):
    """sum_{n>=1} exp(-pi n^2 x), via Jacobi theta with the modular transform for x < 1"""
    x = mp.mpf(x)
    if x < 1:
        return (mp.jtheta(3, 0, mp.exp(-mp.pi/x))/mp.sqrt(x) - 1)/2
    return (mp.jtheta(3, 0, mp.exp(-mp.pi*x)) - 1)/2

MODELS = {n: Model(n) for n in ('ZETA', 'EPS', 'DH')}

_CACHE = {}
def a_prof(M, w):
    # Phi_M is even; laws17 evaluates the theta sums reliably for u >= 0 only (EPS truncation), so use |u|.
    key = (M.name, mp.nstr(w, 30))
    if key not in _CACHE:
        _CACHE[key] = 2*M.phi(abs(2*mp.log(w)))/w
    return _CACHE[key]

def nodes(name, K=48):
    U = U_MAX[name]
    return [mp.e**(mp.mpf(-U/2) + mp.mpf(U)*j/K) for j in range(K+1)]

def mellin_omega(M, Delta, K=48):
    Delta = mp.mpc(Delta)
    return mp.quad(lambda w: w**(Delta-1)*a_prof(M, w), nodes(M.name, K))

# ---------------- A. the identity A_M(Delta) = xi_M(Delta/2), every Delta in C ----------------
DELTAS = [1, mp.mpc(1, 20), mp.mpc(1, 2*14.134725), mp.mpc(0.4, 7), mp.mpc(1.6, -11),
          mp.mpc(3, 2), mp.mpc(-0.5, 4), 2, 0]
A = {}
for name, M in MODELS.items():
    rows = []
    for D in DELTAS:
        lhs = mellin_omega(M, D)
        rhs = xi_closed(name, mp.mpc(D)/2)
        rows.append(dict(Delta=str(mp.nstr(mp.mpc(D), 10)), A_omega_quad=mp.nstr(lhs, 15),
                         xi_closed_at_Delta_over_2=mp.nstr(rhs, 15),
                         abs_err=mp.nstr(abs(lhs-rhs), 3)))
    A[name] = rows
    print('A', name, max(float(r['abs_err']) for r in rows), flush=True)
OUT['A_identity'] = A

# ---------------- B. KG norm = PSS Plancherel on Re Delta = 1 ----------------
B = {}
for name, M in MODELS.items():
    kg = mp.quad(lambda w: abs(a_prof(M, w))**2*w, nodes(name))
    L = 400 if name == 'ZETA' else 200
    pts = [mp.mpf(-L) + mp.mpf(2*L)*j/80 for j in range(81)]
    pl = mp.quad(lambda nu: abs(xi_closed(name, mp.mpc(1, nu)/2))**2, pts)/(2*mp.pi)
    B[name] = dict(KG_norm=mp.nstr(kg, 15), plancherel_ReDelta1=mp.nstr(pl, 15), rel_err=mp.nstr(abs(kg-pl)/kg, 3))
    print('B', name, B[name], flush=True)
OUT['B_KG_vs_Plancherel'] = B

# ---------------- C. why the factor 2: power maps y = omega^k ----------------
# Lorentz-invariant measure d^3p/2omega = (1/2) omega d omega d^2z (d = 2 celestial dims).
# With y = omega^k:  omega d omega = (1/k) y^{2/k - 1} dy;  Lebesgue in y  <=>  k = 2.
# Mellin: omega^{Delta-1} d omega = (1/k) y^{Delta/k - 1} dy  =>  s = Delta/k.
C = dict(rule="y = omega^k gives s = Delta/k; KG measure omega d omega = (1/k) y^(2/k-1) dy; "
              "flat (Plancherel line Re s = 1/2) iff k = 2, i.e. Delta = 2s; then shadow Delta->2-Delta = s->1-s.")
# k = 1 comparison: integer-frequency tower (Planck / linear theta) gives zeta(Delta) with centre Re Delta = 1/2
pl_rows = []
for D in [mp.mpc(3, 1), mp.mpc(1.5, 4), mp.mpc(2.2, -3)]:
    lhs = mp.quad(lambda w: w**(D-1)/mp.expm1(w), [0, 1, 5, 20, 80])
    pl_rows.append(dict(Delta=mp.nstr(D, 6), planck_mellin=mp.nstr(lhs, 14), Gamma_zeta=mp.nstr(mp.gamma(D)*mp.zeta(D), 14)))
th_rows = []
for D in [mp.mpc(3, 1), mp.mpc(1.5, 4)]:
    lhs = mp.quad(lambda w: w**(D-1)*theta_tail(w*w), [0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.6, 1, 2, 3, 5, 8])
    th_rows.append(dict(Delta=mp.nstr(D, 6), linear_theta_mellin=mp.nstr(lhs, 14), half_GammaR_zeta=mp.nstr(mp.pi**(-D/2)*mp.gamma(D/2)*mp.zeta(D)/2, 14)))
C['k1_planck'] = pl_rows; C['k1_linear_theta'] = th_rows
C['k1_verdict'] = "k=1: Gamma(Delta) zeta(Delta) / Gamma_R(Delta) zeta(Delta): FE Delta <-> 1-Delta, critical line Re Delta = 1/2, NOT the d=2 principal series Re Delta = 1 (PS eq. 4.7: divergent unless c = d/2 = 1)."
print('C', pl_rows[0], th_rows[0], flush=True)
OUT['C_power_maps'] = C

# ---------------- D. raw lattice graviton: soft poles ----------------
def r_zeta(w):   # sum_{n>=1} exp(-pi n^2 omega^4)
    return theta_tail(w**4)
D_ = {}
rows = []
for D in [mp.mpc(3, 2), mp.mpf(5), mp.mpc(2.5, -6)]:
    lhs = mp.quad(lambda w: w**(D-1)*r_zeta(w), [0, 0.3, 0.6, 1, 1.5, 3])
    rhs = mp.pi**(-D/4)*mp.gamma(D/4)*mp.zeta(D/2)/4
    rows.append(dict(Delta=mp.nstr(D, 6), mellin=mp.nstr(lhs, 14), quarter_GammaR_zeta=mp.nstr(rhs, 14)))
D_['ZETA_raw_identity'] = rows
D_['ZETA_raw_small_omega'] = [dict(omega=w, r_minus_soft=mp.nstr(r_zeta(mp.mpf(w)) - (1/(2*mp.mpf(w)**2) - mp.mpf(1)/2), 5))
                              for w in (0.4, 0.5, 0.6, 0.7)]
# residues of (1/4) Gamma_R(Delta/2) zeta(Delta/2) at Delta = 2, 0 and regularity at 1, -1, -2
F = lambda D: mp.pi**(-D/4)*mp.gamma(D/4)*mp.zeta(D/2)/4
eps = mp.mpf('1e-10')
D_['ZETA_residues'] = {str(D0): mp.nstr(eps*F(D0+eps), 10) for D0 in (2, 0)}
D_['ZETA_regular_at'] = {str(D0): mp.nstr(F(mp.mpf(D0)+mp.mpf('1e-15')), 10) for D0 in (1, -1, -2, -3)}
# EPS raw profile: sum' r_Q exp(-pi Q omega^2 / sqrt5)  ~ 1/omega^2 - 1
ME = MODELS['EPS']
def r_eps(w):
    """sum'_{(m,n)} exp(-pi (m^2+5n^2) omega^2/sqrt5) = theta(x/sqrt5) theta(sqrt5 x) - 1,  x = omega^2"""
    x = w*w; s5 = mp.sqrt(5)
    return (1 + 2*theta_tail(x/s5))*(1 + 2*theta_tail(s5*x)) - 1
D_['EPS_raw_small_omega'] = [dict(omega=w, r_minus_soft=mp.nstr(r_eps(mp.mpf(w)) - (1/mp.mpf(w)**2 - 1), 5)) for w in (0.4, 0.6, 0.8)]
D_['EPS_raw_identity'] = []
for D in [mp.mpc(3, 2), mp.mpf(5)]:
    lhs = mp.quad(lambda w: w**(D-1)*r_eps(w), [0, 0.2, 0.5, 1, 2, 4, 8])
    s = D/2
    rhs = (mp.sqrt(5)/mp.pi)**s*mp.gamma(s)*(mp.zeta(s)*mp.dirichlet(s, CHI_M20) + mp.dirichlet(s, CHI_M4)*mp.dirichlet(s, CHI_5))/2
    D_['EPS_raw_identity'].append(dict(Delta=mp.nstr(D, 6), mellin=mp.nstr(lhs, 14), closed=mp.nstr(rhs, 14)))
print('D', D_['ZETA_residues'], D_['EPS_raw_identity'][0], flush=True)
OUT['D_soft_poles'] = D_

# ---------------- E. special graviton points under Delta = 2s ----------------
E = {}
for name in MODELS:
    E[name] = {'A(1)=xi(1/2) [Goldstone/memory, self-shadow]': mp.nstr(xi_closed(name, mp.mpf(1)/2), 12),
               'A(2)=xi(1) [stress tensor, shadow of Delta=0]': mp.nstr(xi_closed(name, mp.mpf(1)), 12),
               'A(0)=xi(0) [subleading soft / superrotation]': mp.nstr(xi_closed(name, mp.mpf(0)), 12),
               'A(-1)=xi(-1/2) [sub-subleading]': mp.nstr(xi_closed(name, mp.mpf(-1)/2), 12)}
OUT['E_special_points'] = E
print('E', E, flush=True)

# ---------------- F. EPS: an invisible graviton primary off the principal series ----------------
rho = mp.findroot(lambda s: xi_closed('EPS', s), mp.mpc(0.933, 15.67))   # Newton start = Round 17 label only
Dz = 2*rho
Aq = mellin_omega(ME, Dz, K=96)
OUT['F_EPS_ghost'] = dict(rho=mp.nstr(rho, 12), Delta=mp.nstr(Dz, 12), ReDelta=mp.nstr(Dz.real, 8),
                          A_omega_quad_at_Delta=mp.nstr(abs(Aq), 3), A_omega_quad_scale_at_ReDelta1=mp.nstr(abs(mellin_omega(ME, mp.mpc(1, Dz.imag), K=96)), 6),
                          note="graviton primary of weight Delta (Re Delta != 1, not delta-normalizable, PS eq. 4.7) has zero overlap with the normalizable EPS graviton state")
print('F', OUT['F_EPS_ghost'], flush=True)

OUT['meta'] = dict(dps=mp.mp.dps, seconds=round(time.time()-t0, 1))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r17b_graviton.json'), 'w'), indent=1)
print('done', OUT['meta'])
