"""Queen's probe (team/queen dh_sinh.py, 047ba9f) applied to the Epstein fake. NUMERICAL (mpmath quad + float64 coefficients).
Window convention (same as Queen): f supported on [-L, L] (L = HALF-width), g = f * f~ supported on [-2L, 2L], so atoms n < e^{2L}.
Vector: f(x) = e^{-i g0 x} sinh(b x), with b = beta0 = Re(rho0) - 1/2 and g0 = Im(rho0), at a zero rho0 of Z_Q.
Exact Levy-Khintchine window form (no Galerkin):
  Q(f) = Omega(0)||f||^2 + int_0^inf (||f||^2 - Re G(u)) nu(u) du - 2 sum_{n<e^{2L}} c(n) n^{-1/2} Re G(log n) + Pole(f).
For Lambda_Q(s) = (sqrt|D|/2pi)^s Gamma(s) Z_Q(s):
  Omega_E(t) = log(|D|/4pi^2) + 2 Re psi(1/2 + i t),  so Omega_E(0) = log(|D|/4pi^2) + 2 psi(1/2),  nu_E(u) = 1/sinh(u/2)
  [from Re psi(1/2+it) - psi(1/2) = int_0^inf e^{-u/2}(1 - cos tu)/(1 - e^{-u}) du].
  c(n) = Dirichlet coefficients of -Z_Q'/Z_Q (principal form, so Z_Q/r(1) = 1 + ...).
  Poles of Lambda_Q at s = 0, 1, exactly as for zeta: Pole(f) = 2 Re(fhat(i/2) conj fhat(-i/2)).
The zeta comparison uses Queen's data unchanged (Omega(0) = psi(1/4) - log pi, nu = e^{u/2}/sinh u, c = Lambda(n)).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import numpy as np, math, mpmath as mp, json, sys
mp.mp.dps = 20
N = 170000   # > e^{12}
def rQ(Q):
    a,b,c = Q; r = np.zeros(N+1); M = int(math.isqrt(4*N))+3
    ys = np.arange(-M, M+1)
    for x in range(-M, M+1):
        v = a*x*x + b*x*ys + c*ys*ys
        v = v[(v > 0) & (v <= N)]
        np.add.at(r, v, 1)
    return r
def logderiv(aa):
    S = np.zeros(N+1); c = np.zeros(N+1)
    for d in range(2, N+1):
        c[d] = aa[d]*math.log(d) - S[d]
        if c[d] != 0.0 and 2*d <= N:
            m = N//d; S[2*d::d] += c[d]*aa[2:m+1]
    return c
lam = np.zeros(N+1); sv = np.ones(N+1, bool); sv[:2] = False
for p in range(2, N+1):
    if sv[p]:
        sv[2*p::p] = False; q = p
        while q <= N: lam[q] = math.log(p); q *= p

def run(L, g0, b, coef, Om0, nu, pole):
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
    return dict(L=L, norm2=nf, Q_over_norm=float(Q/nf), arch=float((Om0*nf+float(arch_int))/nf), prime=float(prime/nf), pole=P/nf)

CASES = [  # (form, zero) ; zeros from ../epstein (CERTIFIED boxes for D=-20, NUMERICAL for D=-23)
  ((1,0,5), 0.9079933666439008, 89.55982797855877),
  ((1,0,5), 0.9329696974854141, 15.668249531278473),
  ((1,1,6), 0.9532604747946607, 16.290215720390393),
]
Ls = [2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
out = []
cache = {}
for Q, x0, y0 in CASES:
    a,b_,c_ = Q; D = b_*b_-4*a*c_
    if Q not in cache:
        r = rQ(Q); aa = r/r[1]; aa[0] = 0; cache[Q] = logderiv(aa)
    cE = cache[Q]
    beta, g0 = x0-0.5, y0
    OmE = float(mp.log(abs(D)/(4*mp.pi**2)) + 2*mp.digamma(0.5))
    for L in Ls:
        rE = run(L, g0, beta, cE, OmE, lambda u: 1/mp.sinh(u/2), True)
        rZ = run(L, g0, beta, lam, float(mp.digamma(0.25)-mp.log(mp.pi)), lambda u: mp.e**(u/2)/mp.sinh(u), True)
        row = dict(form=Q, rho0=[x0,y0], L=L, Q_Epstein_over_norm=rE['Q_over_norm'], Q_zeta_over_norm=rZ['Q_over_norm'], epstein=rE, zeta=rZ)
        out.append(row)
        print(f"{Q} rho0={x0:.4f}+{y0:.3f}i L={L}: Q_E/|f|^2={rE['Q_over_norm']:+.4g} (arch {rE['arch']:+.3f}, prime {rE['prime']:+.3f}, pole {rE['pole']:+.2e}) | Q_zeta/|f|^2={rZ['Q_over_norm']:+.4g}", flush=True)
json.dump(out, open('epstein_sinh.json','w'), indent=1)
