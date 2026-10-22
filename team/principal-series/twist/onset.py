"""Onset of F_theta failure predicted from zeta's even near-null vector. NUMERICAL float64. RH OPEN."""
import math, json, numpy as np
from twist import forms
N = 24
def fmin(L, th):
    Fe, Fo, _ = forms('F', L, N, th); return min(np.linalg.eigvalsh(Fe)[0], np.linalg.eigvalsh(Fo)[0])
def bisect(th, lo=0.3, hi=0.8):
    for _ in range(18):
        mid = (lo + hi) / 2
        if fmin(mid, th) < 0: hi = mid
        else: lo = mid
    return (lo + hi) / 2
# predicted onset: lambda_e(L) = 2*theta^2 * kappa(L)/2  i.e.  Q_ze(v) = Q_zo(theta*u*v) (first order, tanh ~ theta u)
from numpy.polynomial import legendre as leg
def kappa(L):
    Ze, Zo, basis = forms('ZETA', L, 40, 0)
    ez, vz = np.linalg.eigh(Ze)
    # project u*v onto basis via quadrature
    xg, wg = leg.leggauss(200); u = L/2*xg + L/2; w = L/2*wg
    B = basis(u); v = vz[:, 0] @ B
    c = B @ (w * u * v)
    return ez[0], float(c @ Zo @ c)
res = []
for th in [0.25, 0.1, 0.03, 0.01, 0.003]:
    Lb = bisect(th)
    # predicted: smallest L on grid with lam_e(L) <= th^2 * kappa(L)
    Lp = None
    for L in np.arange(0.30, 0.80, 0.01):
        le, k = kappa(L)
        if le <= th**2 * k: Lp = L; break
    r = {"theta": th, "onset_direct_bisect": Lb, "onset_predicted_first_order": Lp}
    print(json.dumps(r)); res.append(r)
json.dump(res, open('onset_results.json', 'w'), indent=1)
