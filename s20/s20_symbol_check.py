"""S20 [A]: check that the S17 completed Weil form on the first prime window equals a single Fourier symbol.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.  RH OPEN, GRH OPEN.
Claim (elementary, PROVED in note): for real g in C_c((-h,h)), h <= h* = log3/2, real primitive chi (conductor f, parity eps):
   Q_chi(g) = int_R |ghat(xi)|^2 sigma_chi(xi) dxi/2pi,
   sigma_chi(xi) = Re psi((1+2eps)/4 + i xi/2) + log(f/pi) - 2 kappa chi(2) cos(xi log2),  kappa = log2/sqrt2.
NUMERICAL cross-check against the S17 Arb matrices (Legendre*(1-u^2)^4 basis, M=12), random vectors; double precision."""
import sys, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "s17")); sys.path.insert(0, os.path.join(HERE, "..", "s14"))
from s17_build import build, forms, arb
from s14_exact import p as basis_poly
import numpy as np, mpmath as mp
M = 12
h = arb(3).log()/2
B = build(M, h, prec=256); Fm = forms(B, h)
hm = float(np.log(3)/2); t = float(np.log(2)); kap = t/np.sqrt(2)
u, w = np.polynomial.legendre.leggauss(400)
P = np.array([np.polyval([float(int(c.p))/int(c.q) for c in basis_poly(j).coeffs()][::-1], u) for j in range(M)])
xi = np.linspace(0, 600, 600001); dx = xi[1]-xi[0]
E = np.exp(-1j*np.outer(xi, hm*u))*w   # (nxi, nu)
tab = np.concatenate([np.linspace(0, 40, 8001), np.linspace(40.05, 600, 11200)])
psi = {e: np.array([float(mp.re(mp.digamma(mp.mpf(1+2*e)/4 + 1j*x/2))) for x in tab]) for e in (0, 1)}
psif = {e: np.interp(xi, tab, psi[e]) for e in (0, 1)}   # fine table + linear interp (check only)
def sigma(f, eps, chi2): return psif[eps] + np.log(f/np.pi) - 2*kap*chi2*np.cos(xi*t)
random.seed(20)
for nm, f, eps, chi2 in [("chi3 (mod 3, odd)", 3, 1, -1), ("chi5 quadratic (mod 5, even)", 5, 0, -1), ("chi4 (mod 4, odd)", 4, 1, 0)]:
    Q = Fm[nm]; N = B["N"]
    Qf = np.array([[float(Q[a, b].mid()) for b in range(M)] for a in range(M)]); Nf = np.array([[float(N[a, b].mid()) for b in range(M)] for a in range(M)])
    for trial in range(2):
        c = np.array([random.uniform(-1, 1) for _ in range(M)])
        G2 = np.abs(hm*(E @ (c @ P)))**2
        wts = np.full(xi.size, dx); wts[0] = wts[-1] = dx/2
        qs = 2*np.sum(wts*G2*sigma(f, eps, chi2))/(2*np.pi); ns = 2*np.sum(wts*G2)/(2*np.pi)
        print(f"{nm:30s} trial {trial}: Q(matrix)={c@Qf@c:.12f}  Q(symbol)={qs:.12f}  diff={abs(c@Qf@c-qs):.1e} ; N: {c@Nf@c:.12f} vs {ns:.12f}")
