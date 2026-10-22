"""Queen S1b: Lemma 5.2's own test vectors f_t = phi - tau_t phi (phi smooth bump on (0, ell)), for DH and zeta.
Q(f_t) = 2C(0) - 2C(t), C(t) = W(tau_t g), g = phi * phi~ (supp [-ell, ell]); f_t lives in a window of half-width L = (ell + t)/2.
For t > ell: C(t) = arch(tau_t g) - sum_n c(n) n^{-1/2} g(log n - t);  arch(tau_t g) = -int g(x - t) nu(|x|) dx  (Levy-Khintchine kernel off 0).
C(0) = Omega(0)||phi||^2 + int_0^inf (||phi||^2 - g(u)) nu(u) du (+ pole for zeta: 2 phihat(i/2) phihat(-i/2); no prime atoms when ell < log 2).
Rayleigh R(t) = Q(f_t)/||f_t||^2 = (C(0) - C(t))/||phi||^2 for t > ell. NUMERICAL (float64, trapezoid quadrature)."""
import numpy as np, math, json
N = 400000
xi = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
per = {0: 0.0, 1: 1.0, 2: xi, 3: -xi, 4: -1.0}
a = np.array([per[n % 5] for n in range(N+1)]); a[0] = 0
S = np.zeros(N+1); c = np.zeros(N+1)
for d in range(2, N+1):
    c[d] = a[d]*math.log(d) - S[d]
    if c[d] != 0.0 and 2*d <= N:
        m = N//d; S[2*d::d] += c[d]*a[2:m+1]
lam = np.zeros(N+1); sv = np.ones(N+1, bool); sv[:2] = False
for p in range(2, N+1):
    if sv[p]:
        sv[2*p::p] = False; q = p
        while q <= N: lam[q] = math.log(p); q *= p
ell = 0.5; h = 1e-4
x = np.arange(h/2, ell, h); phi = np.exp(-1/(x*(ell - x)))
phi /= math.sqrt(np.sum(phi**2)*h)          # ||phi|| = 1
gfull = np.convolve(phi, phi[::-1])*h        # lags -(ell-h) .. (ell-h)
lags = (np.arange(len(gfull)) - (len(phi) - 1))*h
gfun = lambda u: np.interp(u, lags, gfull, left=0.0, right=0.0)
import mpmath as _mp
def quad(f,a,b,limit=None):
    return (float(_mp.quad(lambda u: float(f(float(u))), [a,b] if b!=np.inf else [a, a+5, a+40, _mp.inf])),)
def C0(Om0, nu, pole):
    I = quad(lambda u: (1 - gfun(u))*nu(u), 0, ell, limit=400)[0] + quad(nu, ell, np.inf)[0]
    P = 0.0
    if pole:
        P = 2*np.sum(phi*np.exp(-x/2))*h*np.sum(phi*np.exp(x/2))*h
    return Om0 + I + P
n = np.arange(2, N+1); ln = np.log(n)
def Ct(coef, t, nu):
    m = np.abs(ln - t) < ell
    prime = np.sum(coef[2:][m]*n[m]**-0.5*gfun(ln[m] - t))
    uu = np.linspace(t - ell, t + ell, 401); arch = -np.trapezoid(gfun(uu - t)*nu(uu), uu)
    return arch - prime
nuZ = lambda u: np.exp(u/2)/np.sinh(u); nuD = lambda u: np.exp(-u/2)/np.sinh(u)
import mpmath as mp
c0Z = C0(float(mp.digamma(0.25) - mp.log(mp.pi)), nuZ, True); c0D = C0(float(mp.digamma(0.75) - mp.log(mp.pi/5)), nuD, False)
ts = np.arange(1.0, math.log(N) - ell, 0.002)
CZ = np.array([Ct(lam, t, nuZ) for t in ts]); CD = np.array([Ct(c, t, nuD) for t in ts])
RZ = c0Z - CZ; RD = c0D - CD
out = dict(C0_zeta=c0Z, C0_DH=c0D, rows=[])
for T in np.arange(2.0, 12.6, 0.5):
    m = ts <= T
    out["rows"].append(dict(t_max=float(T), L=float((ell+T)/2), DH_minR=float(RD[m].min()), DH_t=float(ts[m][RD[m].argmin()]),
                            ZETA_minR=float(RZ[m].min()), ZETA_t=float(ts[m][RZ[m].argmin()])))
    print(out["rows"][-1], flush=True)
json.dump(out, open("dh_lemma_vectors.json", "w"), indent=1)
print("C0", c0Z, c0D)
