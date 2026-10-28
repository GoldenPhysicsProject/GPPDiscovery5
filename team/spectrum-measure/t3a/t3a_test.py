# NUMERICAL (float64; evidence, not proof). Script behind the Lambda = 2-5 table in team/spectrum-measure/T3a_plunge_mellin.md (65bba74).
# Usage: python3 t3a_test.py 2 3 4 5   (Lambda values). RH is OPEN; this script is RH-silent (archimedean prolates only).
# Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
# NUMERICAL test of one step of Queen (T3a): low-Mellin fraction of transition-band (plunge) even prolates.
import numpy as np, sys
from numpy.polynomial.legendre import legval, leggauss
def prolates(c, N):
    # commuting operator W = -d(1-x^2)d + c^2 x^2 in normalized Legendre basis P_k sqrt(k+1/2)
    k = np.arange(N); A = np.diag(k*(k+1) + c*c*(2*k*(k+1)-1)/((2*k+3)*(2*k-1)))
    off = c*c*(k[:-2]+1)*(k[:-2]+2)/((2*k[:-2]+3)*np.sqrt((2*k[:-2]+1)*(2*k[:-2]+5)))
    A += np.diag(off, 2) + np.diag(off, -2)
    w, V = np.linalg.eigh(A); return w, V
def evalf(V, n, x):
    N = V.shape[0]; coef = V[:, n]*np.sqrt(np.arange(N)+0.5); return legval(x, coef)
X, Wq = leggauss(4000)
for Lam in [float(a) for a in sys.argv[1:]]:
    c = 2*np.pi*Lam**2; N = int(1.5*c) + 200
    w, V = prolates(c, N); S = np.log(Lam); ss = np.linspace(-S, S, 81)
    xs = (X+1)/2; wx = Wq/2  # [0,1]
    ys = xs; out = []
    for n in range(0, int(2*c/np.pi)+40, 2):  # even states
        f0 = evalf(V, n, xs); norm = np.sum(wx*f0**2)  # int_0^1 psi^2
        # concentration lambda: mu psi(x) = int e^{icxy} psi(y) dy ; at x=xs[k] best with large |psi|
        j = np.argmax(np.abs(f0)); fy = evalf(V, n, X)
        mu = np.sum(Wq*np.exp(1j*c*xs[j]*X)*fy)/f0[j]; lam = c/(2*np.pi)*abs(mu)**2
        g = evalf(V, n, ys**2)
        M = np.array([np.sum(wx*g*2*ys**(2j*s)) for s in ss])
        frac = np.trapezoid(np.abs(M)**2, ss)/(2*np.pi*norm)
        out.append((n, lam, frac, evalf(V,n,np.array([0.0]))[0]**2/norm))
    plunge = [o for o in out if 1e-12 < 1-o[1] and o[1] >= 0.5]
    keep = [o for o in out if 1-o[1] <= 1e-12]
    print("Lambda=%.2f c=%.1f 2c/pi=%.1f S=%.3f" % (Lam, c, 2*c/np.pi, S))
    print("  kept (1-lam<=1e-12): max frac %.3e (n=%d)" % max((o[2], o[0]) for o in keep) if keep else "  none kept")
    for o in plunge: print("  plunge n=%d 1-lam=%.2e lowMellinFrac=%.3e psi(0)^2/norm=%.3f  frac*Lam^2/S=%.3f" % (o[0], 1-o[1], o[2], o[3], o[2]*Lam**2/S))
