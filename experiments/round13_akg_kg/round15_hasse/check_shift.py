"""Independent grid check of B_s (s>0): W(f) for f = g + tau_s h computed by brute-force quadrature of Q(A_f) on a grid
must equal B_0(g,g) + B_0(h,h) + 2 Bsym_s(g,h)."""
import numpy as np, math, json, hasse_core as hc, weil_galerkin as wg
from scipy.integrate import quad
def Qgrid(mo, x, f):
    h = x[1]-x[0]; A = np.correlate(f, f, 'full')*h; t = (np.arange(len(A)) - (len(f)-1))*h
    Ap = A[t >= -1e-12]; tp = t[t >= -1e-12]; A0 = Ap[0]
    w = mo['w']; T = tp[-1]
    tt = tp[1:]; arch = mo['C']*A0 + 2*np.trapezoid(w(tt)*(A0 - Ap[1:]), tt) + 2*A0*quad(lambda u: float(w(u)), T, np.inf)[0]
    # small-t correction: integrand ~ (A0-A(t))/(2t) -> finite; trapezoid from first node ok to O(h)
    kp = hc.kpole(mo); pole = np.trapezoid(A*kp(t), t)
    prime = 0.0
    for n, c in zip(mo['n'], mo['c']):
        if math.log(n) < T: prime += 2*c/math.sqrt(n)*np.interp(math.log(n), tp, Ap)
    return arch + pole - prime
out = []
for name in ['ZETA', 'EPS', 'FTH', 'FTAU']:
    mo = hc.MODELS[name](); l = 0.6; M = hc.Mof(l, 20.0)
    for s in [0.45, math.log(2), math.log(3), 1.7]:
        rng = np.random.default_rng(1); g = rng.normal(size=M); hh = rng.normal(size=M)
        x = np.linspace(-l/2, l/2 + s, 40001)
        phi = lambda y: wg.basis(l, M, y)
        f = g @ phi(x) + hh @ phi(x - s)
        B0 = hc.Bs(mo, l, M, 0.0)['B']; B0 = 0.5*(B0+B0.T); Bsh = hc.Bs(mo, l, M, s)['B']
        pred = g@B0@g + hh@B0@hh + g@Bsh@hh + hh@Bsh.T@g
        got = Qgrid(mo, x, f)
        out.append(dict(model=name, s=s, pred=float(pred), grid=float(got), rel=float(abs(pred-got)/abs(pred))))
        print(out[-1], flush=True)
json.dump(out, open('check_shift.json', 'w'), indent=1)
