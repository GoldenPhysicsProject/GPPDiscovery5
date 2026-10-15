"""S20c [A]: (UNTRUSTED float LP) test whether the POINTWISE worst-phase symbol
   sigma_low(xi) = Re psi((1+2eps)/4 + i xi/2) + log(f/pi) - 2 kappa   (<= sigma_z for every |z|<=1)
admits a dual cosine certificate. If yes at (f, eps), it covers every primitive chi of parity eps and conductor >= f, complex g."""
import numpy as np, scipy.special as sp, sys
from scipy.optimize import linprog
t = np.log(2); kap = t/np.sqrt(2); L = np.log(3)
for f, eps, U, D, B in [(12, 0, 10, 0.1, 4.5), (7, 1, 10, 0.1, 3.0), (12, 0, 10, 0.1, 6.0), (7, 1, 10, 0.1, 4.5)]:
    K = int(round((U-L)/D))+1; u = L+D*np.arange(K); nph = 8; Bp = B/np.cos(np.pi/nph)
    sl = lambda x: sp.psi((1+2*eps)/4+0.5j*x).real+np.log(f/np.pi)-2*kap
    xs = np.arange(0, 5000, 0.01); Xi = xs[np.argmax(sl(xs) > Bp+0.05)]
    xi = np.arange(0, Xi, 0.01); s = sl(xi)
    A1 = np.hstack([np.ones((xi.size, 1)), np.cos(np.outer(xi, u))])
    th = np.linspace(0, 2*np.pi, 800, endpoint=False)
    A2 = np.hstack([np.zeros((800*nph, 1)), np.vstack([np.cos(np.outer(th, np.arange(K))+ph) for ph in 2*np.pi*np.arange(nph)/nph])])
    r = linprog(np.r_[-1, np.zeros(K)], A_ub=np.vstack([A1, A2]), b_ub=np.r_[s, np.full(800*nph, B)], bounds=[(None, None)]*(K+1), method="highs-ds")
    print(f"f={f} eps={eps} B={B} K={K} Xi={Xi:.0f}: c_LP={-r.fun:.6f}", flush=True)
