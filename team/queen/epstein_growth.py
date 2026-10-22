"""Queen S1b: window Weil form of an Epstein zeta Z_Q(s) = sum' Q(m,n)^{-s} (disc -20) on f = e^{-i g0 x} sinh(b x),
f in L^2[-L,L] (half-width L). Normalise Dt(s) = (n0^s / r(n0)) Z_Q(s) = 1 + E(s), E = sum_{n>n0} (r(n)/r(n0)) (n/n0)^{-s}
(generalised Dirichlet series, frequencies log(n/n0) > 0). FE: (sqrt(disc)/(2 pi n0))^s Gamma(s) Dt(s) symmetric under s<->1-s,
simple poles of the completed function at s = 0, 1.  -Dt'/Dt = sum_lambda c_lambda e^{-lambda s}, with c solving
(1+E)*c = lambda E as measures (grid step h, blockwise causal FFT recursion).
Q = Omega(0)||f||^2 + int (||f||^2 - Re G(u)) nu(u) du - 2 sum c_l e^{-l/2} Re G(l) + pole,
Omega(0) = 2 log(sqrt(disc)/(2 pi n0)) + 2 psi(1/2), nu(u) = 1/sinh(u/2). NUMERICAL."""
import numpy as np, math, mpmath as mp, json
mp.mp.dps = 20
def run_form(a, bq, cq, b, g0, Ls, h=1e-5):
    lam_max = 2*max(Ls)
    n0 = min(a*x*x + bq*x*y + cq*y*y for x in range(-3, 4) for y in range(-3, 4) if x or y)
    N = int(n0*math.exp(lam_max)) + 1
    r = np.zeros(N+1)
    disc = 4*a*cq - bq*bq
    Y = int(math.sqrt(4*a*N/disc)) + 2
    for y in range(-Y, Y+1):
        A = a; B = bq*y; C = cq*y*y
        D = B*B - 4*A*(C - N)
        if D < 0: continue
        x1 = math.floor((-B - math.sqrt(D))/(2*A)) - 1; x2 = math.ceil((-B + math.sqrt(D))/(2*A)) + 1
        xs = np.arange(x1, x2+1); v = A*xs*xs + B*xs + C
        v = v[(v >= 1) & (v <= N)]
        np.add.at(r, v, 1)
    r0 = r[n0]
    M = int(lam_max/h) + 2
    E = np.zeros(M)
    ns = np.nonzero(r)[0]; ns = ns[ns > n0]
    idx = np.rint(np.log(ns/n0)/h).astype(int); keep = idx < M
    np.add.at(E, idx[keep], r[ns[keep]]/r0)
    lamgrid = np.arange(M)*h
    c = np.zeros(M); rhs = lamgrid*E
    B = max(1, int(math.log(ns[0]/n0)/h) - 1)
    for s0 in range(0, M, B):
        s1 = min(M, s0+B)
        if s0 > 0:
            n = 1 << int(math.ceil(math.log2(s1 + s0 + 1)))
            conv = np.fft.irfft(np.fft.rfft(E[:s1], n)*np.fft.rfft(c[:s0], n), n)[:s1]
            c[s0:s1] = rhs[s0:s1] - conv[s0:s1]
        else:
            c[s0:s1] = rhs[s0:s1]
    Om0 = 2*math.log(math.sqrt(disc)/(2*math.pi*n0)) + 2*float(mp.digamma(0.5))
    nu = lambda u: 1/mp.sinh(u/2)
    out = []
    for L in Ls:
        G = lambda u: (math.sinh(b*(2*L-u))/(2*b) - (2*L-u)*math.cosh(b*u)/2) if u < 2*L else 0.0
        nf = G(0.0)
        ReG = lambda u: math.cos(g0*u)*G(u)
        pts = [0] + [2*L*j/800 for j in range(1, 801)]
        arch = Om0*nf + float(mp.quad(lambda u: (nf - ReG(float(u)))*nu(u), pts) + nf*mp.quad(nu, [2*L, mp.inf]))
        m = lamgrid < 2*L; lg = lamgrid[m]
        Gk = np.cos(g0*lg)*(np.sinh(b*(2*L-lg))/(2*b) - (2*L-lg)*np.cosh(b*lg)/2)
        prime = 2*np.sum(c[m]*np.exp(-lg/2)*Gk)
        fa = mp.quad(lambda x: mp.e**(-1j*g0*x)*mp.sinh(b*x)*mp.e**(-x/2), mp.linspace(-L, L, 60))
        fb = mp.quad(lambda x: mp.e**(-1j*g0*x)*mp.sinh(b*x)*mp.e**(x/2), mp.linspace(-L, L, 60))
        P = float(2*mp.re(fa*mp.conj(fb)))
        Q = arch - prime + P
        out.append(dict(L=L, Q_over_norm=Q/nf, arch=arch/nf, prime=prime/nf, pole=P/nf, norm2=nf))
        print(dict(form=(a, bq, cq), b=b, g0=g0, **out[-1]), flush=True)
    return out
res = {}
Ls = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
res["x2+5y2_validation"] = run_form(1, 0, 5, 0.908-0.5, 89.56, [2.0, 2.5, 3.0, 4.0])
res["2x2+2xy+3y2"] = run_form(2, 2, 3, 2.07941410481-0.5, 7.62886758778, Ls)
json.dump(res, open("epstein_growth.json", "w"), indent=1)
