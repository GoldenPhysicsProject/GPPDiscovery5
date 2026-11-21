"""Round 22, approach 1: the Fejer-synthesized network m_X, built from primes n <= X plus the archimedean place.
rho_X(g0) := sum_rho h(gamma_rho) with h(r) = K(r-g0) + K(r+g0), K(r) = (L/2pi) sinc^2(rL/2), L = log X (Fejer kernel,
K >= 0, int K = 1, Fourier support [-L, L]). By the explicit formula (r22_models docstring) rho_X is a FINITE expression:
  rho_X(g0) = pole terms + archimedean integral - (2/pi) sum_{n<=X} Lam(n) n^{-1/2} (1 - log n/L) cos(g0 log n).
The network: m_X(u) = int_0^inf rho_X(g) dg / (u + g^2)   (Stieltjes iff rho_X >= 0; m_X -> m as X -> inf).
Under RH rho_X = sum over real gamma of K-bumps >= 0 for every X; an off-line zero makes rho_X < 0 somewhere for large X.
No zeros are inputs; ZETA zeros (mpmath) are used ONLY to validate the code (zero side vs prime side)."""
import sys, json, math, numpy as np
from r22_models import MODELS, lam_table
EG = 0.5772156649015329
XS = [10, 30, 100, 300, 1000, 3000, 10000]
GMAX = dict(ZETA=100.0, FTH=60.0, EPS=60.0, DH=100.0)
xg, wg = np.polynomial.legendre.leggauss(12)

def K(w, L):
    y = w*L/2
    return (L/(2*np.pi))*np.where(np.abs(y) < 1e-12, 1.0, (np.sin(y)/np.where(y == 0, 1, y))**2)

def arch_I(a, b, L, g0):
    """(1/2pi) int h(r) psi(a + i b r) dr = -EG g(0) + int_0^inf (e^{-t} g(0) - e^{-a t} g(b t))/(1-e^{-t}) dt"""
    g00 = 1/np.pi; T = L/b
    edges = np.linspace(0, T, int(T/0.02) + 2)
    t = ((edges[:-1, None] + edges[1:, None])/2 + (edges[1:, None] - edges[:-1, None])/2*xg[None, :]).ravel()
    w = (((edges[1:] - edges[:-1])/2)[:, None]*wg[None, :]).ravel()
    gbt = np.cos(np.outer(g0, b*t))*(1 - b*t/L)/np.pi                 # (G, T)
    integ = (np.exp(-t)*g00 - np.exp(-a*t)*gbt)/(-np.expm1(-t))
    return -EG*g00 + integ@w - g00*math.log(-math.expm1(-T))

def rho(name, X, g0, lam):
    M = MODELS[name]; L = math.log(X)
    out = 2*(1/np.pi)*M['logQ']*np.ones_like(g0)
    for (lmb, a, b) in M['arch']:
        out += 2*lmb*arch_I(a, b, L, g0)
    for z in M['poles']:
        r = -1j*(z - 0.5)
        out += 2*(K(r - g0, L) + K(r + g0, L)).real
    n = np.arange(2, int(X) + 1); lm = lam[2:int(X) + 1]; nz = lm != 0
    n = n[nz]; c = lm[nz]*n**-0.5*(1 - np.log(n)/L)
    out -= (2/np.pi)*(np.cos(np.outer(g0, np.log(n)))@c)
    return out

if __name__ == '__main__':
    models = sys.argv[1].split(',') if len(sys.argv) > 1 else ['ZETA', 'FTH', 'EPS', 'DH']
    if len(sys.argv) > 2:      # extended run: X list and gamma range from the command line
        XS[:] = [int(float(x)) for x in sys.argv[2].split(',')]; GMAX[models[0]] = float(sys.argv[3])
    res = {}
    for name in models:
        lam = lam_table(name, max(XS) + 1)
        g0 = np.arange(0.0, GMAX[name] + 1e-9, 0.01)
        row = {}
        for X in XS:
            r = np.concatenate([rho(name, X, g0[i:i + 200], lam) for i in range(0, len(g0), 200)])
            k = int(np.argmin(r))
            neg = g0[r < 0]
            row[str(X)] = dict(min=float(r[k]), at=float(g0[k]), frac_negative=float(np.mean(r < 0)),
                               first_negative_gamma=float(neg[0]) if len(neg) else None, mean=float(r.mean()))
            print(name, X, row[str(X)], flush=True)
        res[name] = row
        json.dump(res, open(f'r22_fejer_{"_".join(models)}{"_ext" if len(sys.argv) > 2 else ""}.json', 'w'), indent=1)
