"""Round 21 shared models: the d = 1 integer states F_M in L^2(R_{>0}, dx) and log-coordinate kernels.
x = e^u, (Vf)(u) = e^{u/2} f(e^u) is unitary L^2(dx) -> L^2(du); the shadow J f(x) = x^{-1} f(1/x) becomes u -> -u.
phi_M(u) = e^{u/2} F_M(e^u) is even (J F_M = F_M) and  int phi_M(u) e^{(s-1/2)u} du = M_M(s)  (strip).
Built from integers / Dirichlet coefficients / Gamma factors only.  No zeta zeros are used."""
import math, numpy as np, mpmath as mp

# ---------------- numpy (double) versions, used for Hankel spectra ----------------
def _thplus(y, nmax=40):            # sum_{n>=1} e^{-pi n^2 y}, y > 0 (use inversion for small y)
    y = np.asarray(y, float); out = np.zeros_like(y)
    big = y >= 1.0; sm = ~big
    n = np.arange(1, nmax + 1)[:, None]
    if big.any(): out[big] = np.exp(-np.pi * n**2 * y[big][None, :]).sum(0)
    if sm.any():
        ys = y[sm]; th = 1 + 2 * np.exp(-np.pi * n**2 / ys[None, :]).sum(0)
        out[sm] = (th / np.sqrt(ys) - 1) / 2
    return out

def F_zeta(x):  x = np.asarray(x, float); return _thplus(x * x) - 1 / (2 * x)
def F_eps(t):   # x^2 + 5 y^2, Mellin = (sqrt5/pi)^s Gamma(s) Z_Q(s)
    t = np.asarray(t, float); r5 = math.sqrt(5)
    return (1 + 2 * _thplus(t / r5)) * (1 + 2 * _thplus(r5 * t)) - 1 - 1 / t
KAP = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
DH_A = {0: 0.0, 1: 1.0, 2: KAP, 3: -KAP, 4: -1.0}          # Davenport-Heilbronn coefficients (periodic mod 5)
def F_dh(x, nmax=200):  # odd-type theta: Mellin = 1/2 (5/pi)^{(s+1)/2} Gamma((s+1)/2) L_DH(s)
    x = np.asarray(x, float); sh = x.shape; x = np.minimum(x.ravel(), 1e6)
    n = np.arange(1, nmax + 1); a = np.array([DH_A[k % 5] for k in n])
    return (a[:, None] * n[:, None] * x[None, :] * np.exp(-np.pi * n[:, None]**2 * x[None, :]**2 / 5)).sum(0).reshape(sh)

def phi_from_F(F):
    return lambda u: np.exp(np.asarray(u) / 2) * F(np.exp(np.asarray(u, float)))

def phi_even(F):
    p = phi_from_F(F); return lambda u: p(np.abs(np.asarray(u, float)))

def phi_fth(theta=0.2, V=110.0, h=0.004):
    """F_theta state: M(s) = M_zeta(s+theta) M_zeta(s-theta) (Gamma_R zeta(s-theta) zeta(s+theta)/4).
    phi_theta = (e^{theta u} phi_zeta) * (e^{-theta u} phi_zeta)  (additive convolution).  Tabulated on a grid."""
    pz = phi_even(F_zeta)
    v = np.arange(-V, V + h / 2, h); a = pz(v) * np.exp(theta * v); b = pz(v) * np.exp(-theta * v)
    from scipy.signal import fftconvolve
    c = fftconvolve(a, b) * h                     # grid offset: index k <-> u = -2V + k h
    u = -2 * V + np.arange(c.size) * h
    m = (u >= -0.0001) & (u <= V)
    uu, cc = u[m], c[m]
    return lambda w: np.interp(np.abs(np.asarray(w, float)), uu, cc)

# ---------------- mpmath versions, used for transmission / boundary-flux checks ----------------
CHI_M4 = (0, 1, 0, -1); CHI_5 = (0, 1, -1, -1, 1)
CHI_M20 = [CHI_M4[n % 4] * CHI_5[n % 5] for n in range(20)]
def ZQ(s): return mp.zeta(s) * mp.dirichlet(s, CHI_M20) + mp.dirichlet(s, CHI_M4) * mp.dirichlet(s, CHI_5)
def GR(s): return mp.pi**(-s / 2) * mp.gamma(s / 2)
def M_zeta(s): return GR(s) * mp.zeta(s) / 2
def M_eps(s):  return (mp.sqrt(5) / mp.pi)**s * mp.gamma(s) * ZQ(s)
def mp_thplus(y):
    y = mp.mpf(y)
    if y < 1: return (mp.jtheta(3, 0, mp.exp(-mp.pi / y)) / mp.sqrt(y) - 1) / 2
    return (mp.jtheta(3, 0, mp.exp(-mp.pi * y)) - 1) / 2
def mpF_zeta(x): return mp_thplus(x * x) - 1 / (2 * x)
def mpF_eps(t):
    r5 = mp.sqrt(5); return (1 + 2 * mp_thplus(t / r5)) * (1 + 2 * mp_thplus(r5 * t)) - 1 - 1 / t
