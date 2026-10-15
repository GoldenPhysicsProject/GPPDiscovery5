"""S20 [G] NUMERICAL: chi7 cubic symbol Ritz test (S17 basis, M=12), real-g (even part of sigma) vs complex-g (full sigma).
Real-g value should be consistent with S19b/Astra (even-g Ritz 0.304 at M=16).  Double precision; sanity only."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "s14"))
from s14_exact import p as basis_poly
import numpy as np, scipy.special as sp, scipy.linalg as la
M = 12; hm = np.log(3)/2; t = np.log(2); kap = t/np.sqrt(2)
u, w = np.polynomial.legendre.leggauss(400)
P = np.array([np.polyval([float(int(c.p))/int(c.q) for c in basis_poly(j).coeffs()][::-1], u) for j in range(M)])
xi = np.linspace(-400, 400, 160001); dx = xi[1]-xi[0]
G = np.vstack([hm*(np.exp(-1j*np.outer(xc, hm*u))*w) @ P.T for xc in np.array_split(xi, 80)])      # (nxi, M)
ps = sp.psi(0.25+0.5j*xi).real + np.log(7/np.pi)
chi2 = np.exp(2j*np.pi/3)
full = ps - 2*kap*np.real(chi2*np.exp(1j*xi*t)); even = ps + kap*np.cos(xi*t)
def mat(s): return (G.conj().T*(s*dx)) @ G/(2*np.pi)
Nm = mat(np.ones_like(xi))
for nm, s in [("chi7 real g (even part)", even), ("chi7 complex g (full symbol)", full), ("chi7bar complex g", full[::-1])]:
    Q = mat(s); print(f"{nm:30s} Ritz lambda_min = {la.eigh(Q, Nm, eigvals_only=True)[0]:.6f}")
