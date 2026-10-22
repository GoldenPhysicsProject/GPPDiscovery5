"""Queen S1: Q_DH on the off-line-pair vector f_L(x) = e^{-i g0 x} sinh(b x) on [-L,L] (b = 0.308517, g0 = 85.699348).
Exact Levy-Khintchine form of the window Weil form (no Galerkin):
  Q(f) = Omega(0)||f||^2 + int_0^inf (||f||^2 - Re G(u)) nu(u) du - 2 sum_{n<e^{2L}} c(n) n^{-1/2} Re G(log n)  [+ pole; DH has none]
  G(u) = <f, tau_u f> = e^{-i g0 u}[sinh(b(2L-u))/(2b) - (2L-u)cosh(bu)/2] (0<=u<=2L), 0 beyond.
DH Gamma factor Gamma((s+1)/2), conductor 5: Omega_DH(0) = psi(3/4) - log(pi/5), nu(u) = e^{-u/2}/sinh u. Also zeta with the same vector
(Omega(0) = psi(1/4) - log pi, nu = e^{u/2}/sinh u, pole 2Re(fhat(i/2) conj fhat(-i/2))). NUMERICAL (mpmath quad, float coefficients)."""
import numpy as np, math, mpmath as mp, json
mp.mp.dps = 20
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
g0, b = 85.699348, 0.308517
def run(L, coef, Om0, nu, pole):
    G = lambda u: (math.sinh(b*(2*L-u))/(2*b) - (2*L-u)*math.cosh(b*u)/2) if u < 2*L else 0.0
    nf = G(0.0)
    ReG = lambda u: math.cos(g0*u)*G(u)
    pts = [0] + [2*L*j/4000 for j in range(1, 4001)]
    arch_int = mp.quad(lambda u: (nf - ReG(float(u)))*nu(u), pts) + nf*mp.quad(nu, [2*L, mp.inf])
    X = min(int(math.exp(2*L)), N); k = np.arange(2, X+1); lk = np.log(k)
    Gk = np.cos(g0*lk)*(np.sinh(b*(2*L-lk))/(2*b) - (2*L-lk)*np.cosh(b*lk)/2)
    prime = 2*np.sum(coef[2:X+1]*k**-0.5*Gk)
    P = 0.0
    if pole:
        fa = mp.quad(lambda x: mp.e**(-1j*g0*x)*mp.sinh(b*x)*mp.e**(-x/2), mp.linspace(-L, L, 400))
        fb = mp.quad(lambda x: mp.e**(-1j*g0*x)*mp.sinh(b*x)*mp.e**(x/2), mp.linspace(-L, L, 400))
        P = float(2*mp.re(fa*mp.conj(fb)))
    Q = Om0*nf + float(arch_int) - prime + P
    return dict(L=L, norm2=nf, Q_over_norm=Q/nf, arch=(Om0*nf+float(arch_int))/nf, prime=prime/nf, pole=P/nf)
out = {"DH": [], "ZETA": []}
for L in [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]:
    r = run(L, c, float(mp.digamma(0.75) - mp.log(mp.pi/5)), lambda u: mp.e**(-u/2)/mp.sinh(u), False); out["DH"].append(r); print("DH", r, flush=True)
    r = run(L, lam, float(mp.digamma(0.25) - mp.log(mp.pi)), lambda u: mp.e**(u/2)/mp.sinh(u), True); out["ZETA"].append(r); print("ZETA", r, flush=True)
json.dump(out, open("dh_sinh.json", "w"), indent=1)
