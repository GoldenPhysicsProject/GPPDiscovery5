"""Grok S91: even prolate spheroidal functions PS_{2m,0}(c, eta), c = 2*pi*mu, mu = lambda^2 (Connes-Consani 2106.01715 sec.3),
via the normalised even-Legendre expansion of the prolate operator -d/de (1-e^2) d/de + c^2 e^2 (float64; Bouwkamp/Slepian method).
Fourier sign chi_m: F(psi_m)(y) = chi_m psi_m(y) on [-lambda, lambda], psi_m(x) = PS_{2m,0}(c, x/lambda), self-dual F (e^{-2 pi i x y})."""
import numpy as np
from numpy.polynomial import legendre as Lg
from scipy.linalg import eigh_tridiagonal
def prolates(mu, M, K=220):
    c = 2 * np.pi * mu; k = np.arange(0, 2 * K, 2).astype(float)
    d = k * (k + 1) + c**2 * (2 * k**2 + 2 * k - 1) / ((2 * k - 1) * (2 * k + 3))
    kk = k[:-1]; e = c**2 * (kk + 1) * (kk + 2) / ((2 * kk + 3) * np.sqrt((2 * kk + 1) * (2 * kk + 5)))
    w, V = eigh_tridiagonal(d, e, select='i', select_range=(0, M - 1))
    coef = np.zeros((M, 2 * K))
    for m in range(M):
        a = V[:, m] * np.sqrt(k + 0.5)          # coefficients on P_k (unnormalised)
        if a[np.argmax(np.abs(a))] < 0: a = -a
        coef[m, ::2] = a
    return c, coef
def ps_eval(coef_m, eta):
    return Lg.legval(eta, coef_m)
def chi(mu, coef_m, ngrid=7):
    lam = np.sqrt(mu); c = 2 * np.pi * mu
    x, w = Lg.leggauss(400)
    om = np.linspace(0.05, 0.6, ngrid)
    F = np.array([lam * np.sum(w * ps_eval(coef_m, x) * np.cos(c * x * o)) for o in om])
    P = ps_eval(coef_m, om)
    return float(F @ P / (P @ P))
if __name__ == "__main__":
    for mu in (5.5, 7.5):
        c, co = prolates(mu, 18)
        print(mu, [round(chi(mu, co[m]), 10) for m in range(16)])
