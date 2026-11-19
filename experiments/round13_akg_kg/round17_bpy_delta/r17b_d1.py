"""Round 17B, d = 1 primary frame (Daniel 2026-11-19): s = Delta literally, H = L^2(R_{>0}, dx),
primaries psi_Delta(x) = x^{-Delta}, delta-normalizable iff Re Delta = 1/2, shadow (Jf)(x) = x^{-1} f(1/x): Delta <-> 1-Delta.

Candidates for "each zero rho yields a normalizable 1D state", built from integer data only:
 (a) BPY state f_M (Round 17):            (psi_Delta, f_M)  = xi_M(Delta)                 all Delta
 (b) Tate local vector g(x) = e^{-pi x^2}: int x^{s-1} g dx = Gamma_R(s)/2                   (no zeros: local factor)
 (c) integer lattice ("plane waves of Z") with the Poisson pole removed:
       ZETA: F(x) = sum_{n>=1} e^{-pi n^2 x^2} - 1/(2x);   int x^{s-1} F dx = Gamma_R(s) zeta(s)/2,  0 < Re s < 1
       EPS : F(t) = Theta_Q(t) - 1 - 1/t;                  int t^{s-1} F dt = (sqrt5/pi)^s Gamma(s) Z_Q(s), 0 < Re s < 1
     F is in L^2(dx) (bounded at 0, ~ -c/x at infinity): normalizable, built from the integers, sees every zero.
 (d) Nyman-Beurling one-sided primary  e_rho(x) = x^{conj(rho)-1} 1_{(0,1)}:  ||e_rho||^2 = 1/(2 Re rho - 1) iff Re rho > 1/2.
Checks: ||F||^2 (quadrature) = (1/2pi) int |M[F](1/2+it)|^2 dt (closed forms) -- Mellin-Plancherel on the PS line;
EPS: the overlap of F_EPS with x^{-(1-rho)} vanishes at an off-line rho (Newton from the Round 17 label, illustration only).
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
mp.mp.dps = 25
from laws17 import Model

def theta_tail(x):
    x = mp.mpf(x)
    if x < 1:
        return (mp.jtheta(3, 0, mp.exp(-mp.pi/x))/mp.sqrt(x) - 1)/2
    return (mp.jtheta(3, 0, mp.exp(-mp.pi*x)) - 1)/2

CHI_M4 = [mp.mpf(v) for v in (0, 1, 0, -1)]
CHI_5 = [mp.mpf(v) for v in (0, 1, -1, -1, 1)]
CHI_M20 = [CHI_M4[n % 4]*CHI_5[n % 5] for n in range(20)]
def ZQ(s): return mp.zeta(s)*mp.dirichlet(s, CHI_M20) + mp.dirichlet(s, CHI_M4)*mp.dirichlet(s, CHI_5)
def GR(s): return mp.pi**(-s/2)*mp.gamma(s/2)

F = {'ZETA': lambda x: theta_tail(x*x) - 1/(2*x),
     'EPS': lambda t: (1 + 2*theta_tail(t/mp.sqrt(5)))*(1 + 2*theta_tail(mp.sqrt(5)*t)) - 1 - 1/t}
M = {'ZETA': lambda s: GR(s)*mp.zeta(s)/2,
     'EPS': lambda s: (mp.sqrt(5)/mp.pi)**s*mp.gamma(s)*ZQ(s)}
# exact end pieces: ZETA F = -1/2 + O(e^{-pi/x^2}/x) near 0 and -1/(2x) + O(e^{-pi x^2}) at oo;
#                   EPS  F = -1 + O(e^{-pi/(sqrt5 t)}/t) near 0 and -1/t + O(e^{-pi t/sqrt5}) at oo.
ENDS = {'ZETA': (mp.mpf('0.04'), mp.mpf(8), mp.mpf(-1)/2, mp.mpf(-1)/2), 'EPS': (mp.mpf('0.02'), mp.mpf(45), mp.mpf(-1), mp.mpf(-1))}
def mellin(name, s):
    e, T, c0, cinf = ENDS[name]
    a, b = mp.log(e), mp.log(T)
    mid = mp.quad(lambda u: mp.exp(s*u)*F[name](mp.exp(u)), mp.linspace(a, b, 40))
    return mid + c0*e**s/s + cinf*(-T**(s-1)/(s-1))
PTS = [0, mp.mpf('0.05'), mp.mpf('0.2'), mp.mpf('0.5'), 1, 2, 5, 20, 100, mp.inf]
t0 = time.time(); OUT = {}
c = {}
for name in F:
    rows = []
    for s in [mp.mpc(0.5, 3), mp.mpc(0.3, 10), mp.mpc(0.8, -7), mp.mpf(0.5)]:
        lhs = mellin(name, s)
        rows.append(dict(s=mp.nstr(s, 6), mellin_quad=mp.nstr(lhs, 14), closed=mp.nstr(M[name](s), 14), err=mp.nstr(abs(lhs-M[name](s)), 3)))
    norm = mp.quad(lambda x: F[name](x)**2, PTS)
    L = 120
    pl = 2*mp.quad(lambda t: abs(M[name](mp.mpc(0.5, t)))**2, [mp.mpf(L)*j/60 for j in range(61)])/(2*mp.pi)
    c[name] = dict(strip_identity=rows, L2_norm=mp.nstr(norm, 15), plancherel_Re_half=mp.nstr(pl, 15), rel_err=mp.nstr(abs(norm-pl)/norm, 3),
                   F_at_0plus=mp.nstr(F[name](mp.mpf('1e-3')), 8), x_times_F_at_large=mp.nstr(50*F[name](mp.mpf(50)), 8))
    print(name, c[name]['rel_err'], [r['err'] for r in rows], flush=True)
OUT['c_lattice_state'] = c
rho = mp.findroot(lambda s: ZQ(s), mp.mpc(0.933, 15.67))   # Newton start = Round 17 label only
ov = mellin('EPS', rho)
OUT['c_EPS_ghost'] = dict(
    rho=mp.nstr(rho, 12), invisible_primary_Delta=mp.nstr(1-rho, 12), overlap_quad=mp.nstr(abs(ov), 3),
    scale_same_height_on_line=mp.nstr(abs(M['EPS'](mp.mpc(0.5, rho.imag))), 6),
    d_one_sided_norm_sq=mp.nstr(1/(2*rho.real-1), 10))
print(OUT['c_EPS_ghost'], flush=True)
# (b) Tate local vector: Gamma_R has no zeros
OUT['b_tate'] = dict(check=[dict(s=mp.nstr(s, 6), mellin=mp.nstr(mp.quad(lambda x: x**(s-1)*mp.exp(-mp.pi*x*x), [0, 1, 4, mp.inf]), 12),
                                 half_GammaR=mp.nstr(GR(s)/2, 12)) for s in (mp.mpc(0.5, 14.134725), mp.mpc(0.7, 3))],
                     note="local Tate vector: Mellin = Gamma_R(s)/2, never zero -> sees no zeta zero; zeros only after the global sum over n (candidate c)")
OUT['meta'] = dict(dps=mp.mp.dps, seconds=round(time.time()-t0, 1))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r17b_d1.json'), 'w'), indent=1)
print('done', OUT['meta'])
