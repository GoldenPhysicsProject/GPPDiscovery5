"""Round 20 core: Taylor coefficients of xi_M at 1/2 (zero-free: from the theta law Phi_M of Round 17, laws17.py),
F_M(u) = Xi_M(sqrt u)/Xi_M(0), m_M(u) = F_M'/F_M, its moments c_j (m = sum (-1)^j c_j u^j), the S-fraction
m = c_0/(1 + a_1 u/(1 + a_2 u/(1 + ...))) and the Cauer (RC) ladder. No zeros are used anywhere.
Models: ZETA, EPS (x^2+5y^2), DH, FTH = xi(s-1/4) xi(s+1/4)."""
import os, sys, json, math
import mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'round17_bpy_delta'))
import laws17 as L

STRIP = dict(ZETA=mp.pi/4, EPS=mp.pi/2, DH=mp.pi/4, FTH=mp.pi/4)

def even_moments(name, K, dps):
    """m_{2k} = int Phi_M u^{2k} du, k=0..K, trapezoid with step chosen for dps digits (strip width d: err ~ e^{-2 pi d/h})."""
    mp.mp.dps = dps + 20
    base = 'ZETA' if name == 'FTH' else name
    d = STRIP[name]
    h = 2*mp.pi*d/((dps + 30)*mp.log(10))
    model = L.Model(base)
    us, ph = L.grid(model, h)
    if name != 'FTH':
        return L.moments_even(us, ph, h, K)
    th = mp.mpf(1)/4
    Ma = L.tilted_moments(us, ph, h, 2*K, -th); Mb = L.tilted_moments(us, ph, h, 2*K, th)
    return L.product_moments(Ma, Mb, K)

def F_coeffs(name, K, dps):
    m = even_moments(name, K, dps)
    return [m[k]/mp.factorial(2*k)/m[0] for k in range(K + 1)]   # F(u) = sum f_k u^k, f_0 = 1

def series_inv(a, n):
    b = [mp.mpf(0)]*n; b[0] = 1/a[0]
    for k in range(1, n):
        b[k] = -mp.fsum(a[j]*b[k - j] for j in range(1, min(k, len(a) - 1) + 1))/a[0]
    return b

def m_series(f, n):
    """m = F'/F as power series, n terms"""
    fp = [(k + 1)*f[k + 1] for k in range(len(f) - 1)]
    inv = series_inv(f, n)
    return [mp.fsum(fp[j]*inv[k - j] for j in range(k + 1)) for k in range(n)]

def s_fraction(mu, nmax):
    """mu = power series of m. Returns c0 and a_1..a_nmax with m = c0/(1 + a1 u/(1 + a2 u/(1+...)))."""
    c0 = mu[0]; g = [x/c0 for x in mu]; out = []
    for k in range(nmax):
        if len(g) < 3: break
        ig = series_inv(g, len(g))
        q = ig[1:]                       # (1/g - 1)/u
        a = q[0]
        out.append(a)
        if a == 0: break
        g = [x/a for x in q]
    return c0, out

def ladder(c0, a):
    """Y = 1/m = G1 + 1/(1/(C1 u) + 1/(G2 + 1/(1/(C2 u) + ...))): shunt conductances G_k, series capacitors C_k."""
    G = [1/c0]; C = []
    for k in range(len(a)//2):
        C.append(G[-1]*a[2*k]); G.append(C[-1]/a[2*k + 1])
    return G, C

def m_eval(f, u):
    """m(u) = F'(u)/F(u) from the Taylor coefficients (F entire)."""
    F = mp.fsum(f[k]*u**k for k in range(len(f)))
    Fp = mp.fsum(k*f[k]*u**(k - 1) for k in range(1, len(f)))
    return Fp/F
