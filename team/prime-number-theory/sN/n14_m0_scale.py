# N.14 NUMERICAL scale estimate of the H-norm margin m0^sigma = inf R(xi)/||xi||_H^2 on parity sector sigma,
# xi in C_c^inf(-L0, L0), L0 = 1/2 log 2, R(xi) = (1/2pi) int (1/4+t^2) Omega(t) |xi^|^2, Omega = Re psi(1/4+it/2) - log pi.
# Rayleigh-Ritz on uniform P1 hats (upper estimates of m0; NUMERICAL, no remainder). Norm: H = H^1_{1/4}(R).
import numpy as np
def psi(z, N=12):
    z = np.asarray(z, complex); acc = np.zeros_like(z)
    for k in range(N): acc -= 1/(z + k)
    w = z + N; w2 = 1/(w*w)
    return acc + np.log(w) - 0.5/w - w2*(1/12 - w2*(1/120 - w2*(1/252 - w2*(1/240 - w2/132))))
Om = lambda t: np.real(psi(0.25 + 0.5j*t)) - np.log(np.pi)
L0 = 0.5*np.log(2)
for N in [20, 40, 80, 160]:
    h = 2*L0/N; x = -L0 + h*np.arange(1, N); n = x.size
    Tmax = 400.0/h; dt = 0.01
    t = np.arange(dt/2, Tmax, dt)
    wR = (0.25 + t**2)*Om(t)*h*h*np.sinc(t*h/2/np.pi)**4*dt/np.pi
    d = np.arange(n)*h
    rowR = np.array([np.sum(wR*np.cos(t*dd)) for dd in d])
    # tail t > Tmax: hat^2 = 16 sin^4(th/2)/(h^2 t^4), Omega ~ log(t/2pi)
    # exact means of sin^4(th/2)cos(kth): 3/8, -1/4, 1/16 (k = 0, 1, 2), else 0 -> keeps the row sums consistent
    tl = (16/h**2)*(np.log(Tmax/(2*np.pi)) + 1)/Tmax/np.pi
    for k, m in [(0, 3/8), (1, -1/4), (2, 1/16)]: rowR[k] += tl*m
    GR = rowR[np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])]
    GH = (np.diag(np.full(n, 2/h)) - np.diag(np.full(n-1, 1/h), 1) - np.diag(np.full(n-1, 1/h), -1)) \
         + 0.25*(np.diag(np.full(n, 2*h/3)) + np.diag(np.full(n-1, h/6), 1) + np.diag(np.full(n-1, h/6), -1))
    Rf = np.eye(n)[::-1]
    out = []
    for sgn in (+1, -1):
        P = (np.eye(n) + sgn*Rf)/2; U, s, _ = np.linalg.svd(P); B = U[:, s > 0.5]
        A = B.T @ GR @ B; C = B.T @ GH @ B
        Ci = np.linalg.inv(np.linalg.cholesky(C))
        out.append(np.linalg.eigvalsh(Ci @ A @ Ci.T)[:2])
    print(f"N={N:4d}  even: m0 <= {out[0][0]: .5f} (next {out[0][1]: .4f})   odd: m0 <= {out[1][0]: .5f} (next {out[1][1]: .4f})")
