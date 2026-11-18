"""Round 14F seventh row: unitary twist F_{i tau} = zeta(s-i tau) zeta(s+i tau) (tau=3) vs non-unitary F_theta (theta=0.2),
in the 14B closed-form Galerkin (../round14b_euler/weil_galerkin.py), own explicit-formula data only (no zeros):
  prime  c(n) = 2 cos(tau log n) Lambda(n)       [F_theta: 2 cosh(theta log n) Lambda(n)]
  arch   w(t) = 2 cos(tau t) w_zeta(t),  C = 2 Re psi(1/4 + i tau/2) - 2 log pi
  poles  s = 1 +- i tau -> exponents 1/2 +- i tau: P = 2 int_0^L (A+A^T)(t) * 2cosh(t/2)cos(tau t) dt (Gauss-Legendre)
Pole quadrature validated on zeta (tau=0, kernel cosh(t/2)) against the closed form.
Run from this folder: python3 r14f_twist.py -> r14f_twist.json"""
import sys, os, math, json, numpy as np
from scipy.special import digamma
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.join(HERE, '..', 'round14b_euler'), os.path.join(HERE, '..')]
os.chdir(os.path.join(HERE, '..', 'round14b_euler'))
import weil_galerkin as wg, models14b as md
Z = md.zeta()
XG, WGT = np.polynomial.legendre.leggauss(400)
def pole_quad(L, M, kern):
    P = np.zeros((M, M)); t = 0.5*L*(XG+1); w = 0.5*L*WGT
    for tk, wk in zip(t, w):
        A = wg.Aclosed(L, M, tk); P += wk*kern(tk)*(A+A.T)
    return 2*P
def ftau(tau):
    C = 2*digamma(complex(0.25, tau/2)).real - 2*math.log(math.pi)
    return dict(name=f'FTAU{tau}', n=Z['n'], c=Z['c']*2*np.cos(tau*np.log(Z['n'])), C=C,
                w=lambda t: 2*np.cos(tau*t)*md.wz(t), poles=[])
out = {}
for L in [2.0, 4.0]:
    M = int(math.ceil(45*L/math.pi)) + 8
    a0, p0, q0 = wg.operators(Z, L, M)
    pz = pole_quad(L, M, lambda t: math.cosh(t/2))
    val = float(np.abs(pz - p0).max()/np.abs(p0).max())
    rows = {"M": M, "pole_quadrature_vs_closed_form_on_zeta(rel)": val}
    ez = np.linalg.eigvalsh(a0+p0-q0); rows["zeta"] = {"lam_min": float(ez[0]), "ghosts(<-1e-8*max)": int(np.sum(ez < -1e-8*np.abs(ez).max()))}
    for name, mo, P in [("F_theta 0.2", md.fth(0.2), None), ("F_itau 3", ftau(3.0), "q")]:
        a, p, q = wg.operators(mo, L, M)
        if P == "q": p = pole_quad(L, M, lambda t: 2*math.cosh(t/2)*math.cos(3.0*t))
        ev = np.linalg.eigvalsh(a+p-q)
        rows[name] = {"lam_min": float(ev[0]), "ghosts(<-1e-8*max)": int(np.sum(ev < -1e-8*np.abs(ev).max())), "lam_max": float(ev[-1])}
    out[f"L={L}"] = rows
# twist identity check for theta (digamma closes): psi(1/4+th/2)+psi(1/4-th/2)-2psi(1/4) = -4 int w_z (cosh(th t)-1)
from scipy.integrate import quad
out["twist_identity_theta"] = {str(th): float(digamma(0.25+th/2)+digamma(0.25-th/2)-2*digamma(0.25)
                               + 4*quad(lambda t: md.wz(t)*(math.cosh(th*t)-1), 0, 200, limit=400)[0]) for th in [0.1, 0.2, 0.3]}
out["twist_identity_tau"] = {str(ta): float(2*digamma(complex(0.25, ta/2)).real-2*digamma(0.25)
                             - 4*quad(lambda t: md.wz(t)*(1-math.cos(ta*t)), 0, 200, limit=2000)[0]) for ta in [1.0, 3.0]}
json.dump(out, open(os.path.join(HERE, 'r14f_twist.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
