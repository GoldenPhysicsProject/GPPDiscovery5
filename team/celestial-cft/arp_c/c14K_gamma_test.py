"""C.14.6 numerical test of Conjecture C.14.K (NUMERICAL).
gamma(T, ell) = sup |2P_{>T}(xi_b, xi_s)| / [P_{>T}(xi_b) + P_{>T}(xi_s)] over piecewise-linear xi = xi_b + xi_s,
in the odd and even sectors, = largest canonical correlation (cosine of the minimal angle) between U = {xi_b} and V = {xi_s}
in the P_{>T} geometry, computed via the generalized eigenproblem [[0,G_UV],[G_VU,0]] x = gamma diag(G_UU, G_VV) x.
Compare with gamma* = Omega(T)/(Omega(T) + |Omega(0)|).

Discretization: uniform mesh on [-L, L], step h = 2L/M, strips of n_s elements, so ell = a n_s/(M - n_s) exactly.
U = continuous p.w.-linear functions, arbitrary on the bulk [-L+ell, L-ell], linear (ramp) on each strip, 0 at +-L.
V = p.w.-linear functions vanishing outside the open strips S_- = [-L, -L+ell], S_+ = [L-ell, L] (H^1_0 on each strip).
P(u,v) = (1/2pi) int (1/4+t^2)(Omega(t)-Omega(0)) u^ conj(v^) dt. Position form (exact, from
Omega(t)-Omega(0) = int_R (1 - cos tr) kappa(|r|) dr, kappa(r) = e^{-r/2}/(1-e^{-2r})):
  for p.w.-linear u = sum a_k (x-x_k)_+, v = sum b_m (x-y_m)_+ (compact support),
  P(u,v) = sum a_k b_m K(x_k - y_m),  K(d) = (1/2) Phi1(d) + Psi(d),
  Phi1(d) = int_|d|^inf 2 kappa(r)(r-|d|) dr,
  Psi(d)  = (1/24)[ int_0^|d| kappa(r)(-3|d| r^2) dr + int_|d|^inf kappa(r)(|d|^3 - r^3 - 3 r d^2) dr ].
P_{<T} by Gauss-Legendre in t on [0,T] with hat transform h sinc^2(th/2); P_{>T} = P - P_{<T}.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import numpy as np, mpmath as mp, sys
from scipy.integrate import quad
from scipy.linalg import eigh, toeplitz
a = np.log(2.0)
Om = lambda t: float(mp.re(mp.digamma(0.25 + 0.5j*t)) - mp.log(mp.pi))
Om0 = Om(0)
kap = lambda r: np.exp(-r/2)/(-np.expm1(-2*r))
def Kfun(d):
    d = abs(d)
    if d == 0:
        p1 = quad(lambda r: 2*kap(r)*r, 0, np.inf, limit=400, epsabs=1e-15, epsrel=1e-13)[0]
        ps = quad(lambda r: kap(r)*(-r**3), 0, np.inf, limit=400, epsabs=1e-15, epsrel=1e-13)[0]/24
        return 0.5*p1 + ps
    p1 = quad(lambda r: 2*kap(r)*(r-d), d, np.inf, limit=400, epsabs=1e-15, epsrel=1e-13)[0]
    i1 = quad(lambda r: kap(r)*(-3*d*r*r), 0, d, limit=400, epsabs=1e-15, epsrel=1e-13)[0]
    i2 = quad(lambda r: kap(r)*(d**3 - r**3 - 3*r*d*d), d, np.inf, limit=400, epsabs=1e-15, epsrel=1e-13)[0]
    return 0.5*p1 + (i1+i2)/24
_Kcache = {}
def toeplitz_P(M, h, T, nGL=600):
    # column c_k = P(hat_0, hat_k), k = 0..M
    key = (M, round(h, 15))
    if key not in _Kcache:
        _Kcache[key] = np.array([Kfun(k*h) for k in range(M+3)])
    Kv = _Kcache[key]
    w = np.array([1.0, -2.0, 1.0])/h
    col = np.zeros(M+1)
    for k in range(M+1):
        s = 0.0
        for p in range(3):
            for q in range(3):
                s += w[p]*w[q]*Kv[abs(k + p - q)]
        col[k] = s
    # low band
    x, wt = np.polynomial.legendre.leggauss(nGL)
    t = 0.5*T*(x+1); wt = 0.5*T*wt
    sym = np.array([Om(tt) for tt in t]) - Om0
    hat2 = (h*np.sinc(t*h/(2*np.pi))**2)**2
    base = wt*sym*(0.25 + t*t)*hat2/np.pi
    low = np.array([np.sum(base*np.cos(t*k*h)) for k in range(M+1)])
    return col, col - low
def gamma(M, ns, T, sector, retvec=False):
    L = a*M/(2*(M-ns)) ; ell = 2*L - a; h = 2*L/M
    col, colhi = toeplitz_P(M, h, T)
    n = M - 1                      # interior nodes 1..M-1 (node j at x = -L + j h)
    G = toeplitz(colhi[:n])
    # U basis: bulk nodes j = ns..M-ns (edge nodes ns, M-ns carry ramps)
    EU = []
    for j in range(ns, M-ns+1):
        v = np.zeros(n); v[j-1] = 1.0
        if j == ns:
            for i in range(1, ns): v[i-1] = i/ns
        if j == M-ns:
            for i in range(M-ns+1, M): v[i-1] = (M-i)/ns
        EU.append(v)
    EV = []
    for j in list(range(1, ns)) + list(range(M-ns+1, M)):
        v = np.zeros(n); v[j-1] = 1.0; EV.append(v)
    EU = np.array(EU).T; EV = np.array(EV).T
    R = np.eye(n)[::-1]            # x -> -x
    sgn = 1 if sector == 'even' else -1
    def proj(E):
        F = E + sgn*(R @ E)
        q, r = np.linalg.qr(F)
        keep = np.abs(np.diag(r)) > 1e-10*np.abs(np.diag(r)).max()
        return q[:, keep]
    BU, BV = proj(EU), proj(EV)
    Guu, Gvv, Guv = BU.T@G@BU, BV.T@G@BV, BU.T@G@BV
    # whiten with eigen-decomposition (robust to near-singularity)
    def isqrt(A):
        w, Q = np.linalg.eigh(A); return Q @ np.diag(w**-0.5) @ Q.T, w.min()/w.max()
    Wu, cu = isqrt(Guu); Wv, cv = isqrt(Gvv)
    C = Wu @ Guv @ Wv
    U_, s, Vt = np.linalg.svd(C)
    out = dict(L=L, ell=ell, h=h, gamma=s[0], condU=cu, condV=cv, nU=BU.shape[1], nV=BV.shape[1])
    if retvec:
        xb = BU @ (Wu @ U_[:, 0]); xs = BV @ (Wv @ Vt[0])
        if (xb @ G @ xs) > 0: xs = -xs
        out['xi_b'] = xb; out['xi_s'] = xs
    return out
def gstar(T): return Om(T)/(Om(T) + abs(Om0))
if __name__ == '__main__':
    # sanity: P_full Toeplitz vs frequency quadrature for hat_0 with itself and a neighbour
    M, h = 60, 0.7/60
    col, _ = toeplitz_P(M, h, 14.0)
    for k in (0, 1, 5):
        f = lambda t: (Om(t)-Om0)*(0.25+t*t)*(h*np.sinc(t*h/(2*np.pi))**2)**2*np.cos(t*k*h)/np.pi
        val = float(mp.quad(lambda t: f(float(t)), [0, 50, 200, 1000, 5000, 20000, 1e5]))
        tail = 0  # remaining tail ~ int log t * 16/(h^2 t^4) small
        print(f"sanity k={k}: position {col[k]:.8f}  frequency(0..1e5) {val:.8f}")
    sys.stdout.flush()
