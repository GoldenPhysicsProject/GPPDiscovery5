"""Copied from experiments/round13_akg_kg/round22_synthesis/r22_models.py for Round 25b (fig8).

Round 22 model data (no zeros): local data of each model for the explicit formula.
For xi_M(s) = R(s) Q^s G(s) D(s), D = sum b_n n^{-s}, -D'/D = sum Lam_M(n) n^{-s}, xi_M(s) = xi_M(1-s):
  sum_rho h(gamma_rho) = 2 sum_{z in (1/2, c): R(z)=0} h((z-1/2)/i) + (1/pi) int h(r)[log Q + Re G'/G(1/2+ir)] dr
                         - 2 sum_n Lam_M(n) n^{-1/2} g(log n),   g(x) = (1/2pi) int h(r) e^{irx} dr,
(contour argument over Re s = c > sigma_abs, then FE; R'/R contributes 0 on the line because R(s) = R(1-s)).
G'/G is a sum of terms lam*psi(lam*s + mu); ARCH lists (coeff, a, b) with Re-part integral
  (1/pi) int h(r) Re[lam psi(lam(1/2+ir)+mu)] dr = 2 * lam * (1/2pi) int h(r) psi(a + i b r) dr,  a = lam/2+mu, b = lam."""
import math, numpy as np

TH = 0.25
KAP = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)

MODELS = {
    'ZETA': dict(logQ=-0.5*math.log(math.pi), arch=[(0.5, 0.25, 0.5)], poles=[1.0]),
    'FTH':  dict(logQ=-math.log(math.pi), arch=[(0.5, 0.25 - TH/2, 0.5), (0.5, 0.25 + TH/2, 0.5)], poles=[1 + TH, 1 - TH]),
    'EPS':  dict(logQ=math.log(math.sqrt(20)/(2*math.pi)), arch=[(1.0, 0.5, 1.0)], poles=[1.0]),
    'DH':   dict(logQ=0.5*math.log(5/math.pi), arch=[(0.5, 0.75, 0.5)], poles=[]),
}

def lam_table(name, N):
    """Lam_M(n), n = 0..N (index 0 unused)."""
    n = np.arange(N + 1)
    if name in ('ZETA', 'FTH'):
        lam = np.zeros(N + 1); sieve = np.ones(N + 1, bool); sieve[:2] = False
        for p in range(2, N + 1):
            if sieve[p]:
                sieve[p*p::p] = False
                q = p
                while q <= N: lam[q] = math.log(p); q *= p
        if name == 'FTH':
            lam[1:] = lam[1:]*(n[1:]**TH + n[1:]**(-TH))
        return lam
    if name == 'DH':
        b = np.array([[0, 1, KAP, -KAP, -1][k % 5] for k in range(N + 1)], float)
    else:  # EPS: b_n = r_Q(n)/2, Q = x^2 + 5 y^2
        b = np.zeros(N + 1); M = int(math.isqrt(N)) + 1
        for x in range(-M, M + 1):
            for y in range(-M, M + 1):
                v = x*x + 5*y*y
                if 0 < v <= N: b[v] += 1
        b /= 2
    # Lam = (b log) * b^{-1}  (Dirichlet convolution); b_1 = 1
    binv = np.zeros(N + 1); binv[1] = 1.0; acc = np.zeros(N + 1)
    for d in range(1, N + 1):
        if d > 1: binv[d] = -acc[d]
        if binv[d] != 0 and 2*d <= N:
            acc[2*d::d] += binv[d]*b[2:N//d + 1]
    blog = b*np.log(np.maximum(n, 1))
    lam = np.zeros(N + 1)
    for d in range(1, N + 1):
        if blog[d] == 0: continue
        lam[d::d] += blog[d]*binv[1:N//d + 1]
    return lam
