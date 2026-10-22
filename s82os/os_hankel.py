"""Grok S82-OS: half-window Osterwalder-Schrader (Hankel) form of the Weil distribution vs the full-window Weil form.
NUMERICAL (float64 / numpy). RH OPEN.

Window [-L,L] in log-scale u. Reflection Theta: u -> -u (x -> 1/x). 'Future' test functions h supported in [0,L].
Weil distribution W (S14 normalisation, same as s62/s75 dhritz.py) as a kernel on R:
  ZETA: W(r) = 2cosh(r/2) [pole]  + c0 delta_0 - (1/2)(k1+k2)(|r|) [arch, k1 part regularised at 0]
               - sum_n Lambda(n)/sqrt(n) (delta_{log n} + delta_{-log n})
  DH  : no pole, + log5 delta_0, + (1/2)sech(r/2), prime weights b(n) (S47 comb, kappa from q=5).
For h on [0,L] (orthonormal shifted Legendre basis):
  T(h) = int int h(x)h(y) W(x-y)   (Toeplitz, 'same-side' form)
  B(h) = int int h(x)h(y) W(x+y)   (Hankel = <Theta h, h>_W, the OS form)
Even extension f(u)=h(|u|): Q_W(f)/||f||^2 = (T+B)(h)/||h||^2 ; odd extension: (T-B)(h)/||h||^2.
So Weil positivity on [-L,L] (both parities) <=> -T <= B <= T  <=> T>=0 and C = T^{-1/2} B T^{-1/2} is a contraction.
OS positivity (reflection positivity w.r.t. Theta on E_+ = L^2[0,L]) would be B >= 0.
Usage: os_hankel.py FORM L N   (FORM in ZETA, DH)
"""
import sys, json, math
import numpy as np
from numpy.polynomial import legendre as leg

def mang(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0

def dh_weights(NN):
    s5 = math.sqrt(5); kap = (math.sqrt(10 - 2 * s5) - 2) / (s5 - 1)
    chi = {0: (0, 0), 1: (1, 0), 2: (0, 1), 3: (0, -1), 4: (-1, 0)}
    aa = {n: chi[n % 5][0] + kap * chi[n % 5][1] for n in range(1, NN + 1)}
    bb = {1: 0.0}
    for n in range(2, NN + 1):
        bb[n] = aa[n] * math.log(n) - sum(bb[d] * aa[n // d] for d in range(2, n) if n % d == 0)
    return bb

def build(FORM, L, N, M=96, S=160):
    xg, wg = leg.leggauss(M)
    def basis(x):   # orthonormal shifted Legendre on [0,L]; zero outside
        x = np.asarray(x, float); t = 2 * x / L - 1
        out = np.zeros((N,) + x.shape)
        inside = (x >= 0) & (x <= L)
        for k in range(N):
            c = np.zeros(k + 1); c[k] = 1
            out[k] = np.where(inside, math.sqrt((2 * k + 1) / L) * leg.legval(t, c), 0.0)
        return out
    def conv_plus(s):   # D_ij(s) = int phi_i(x) phi_j(s-x) dx, s in [0,2L]
        a, b = max(0.0, s - L), min(L, s)
        if b <= a: return np.zeros((N, N))
        x = (b - a) / 2 * xg + (a + b) / 2; w = (b - a) / 2 * wg
        P = basis(x); Q = basis(s - x)
        return (P * w) @ Q.T
    def conv_minus(r):  # C_ij(r) = int phi_i(x) phi_j(x-r) dx, r in [0,L]; symmetrised
        a, b = r, L
        if b <= a: return np.zeros((N, N))
        x = (b - a) / 2 * xg + (a + b) / 2; w = (b - a) / 2 * wg
        P = basis(x); Q = basis(x - r)
        Cm = (P * w) @ Q.T
        return 0.5 * (Cm + Cm.T)
    k1 = lambda z: math.exp(-z) / math.sinh(z)
    k2 = lambda z: 2 * math.exp(-z / 4) * math.sinh(3 * z / 4) / math.sinh(z)
    sech = lambda z: 1 / math.cosh(z / 2)
    # e^{+-x/2} moments for the pole
    x = L / 2 * xg + L / 2; w = L / 2 * wg; P = basis(x)
    ap = P @ (w * np.exp(x / 2)); am = P @ (w * np.exp(-x / 2))
    sg, swg = leg.leggauss(S)
    # ---- Hankel B ----
    B = np.zeros((N, N)); comp = {}
    comp['B_pole'] = (np.outer(ap, ap) + np.outer(am, am)) if FORM == 'ZETA' else np.zeros((N, N))
    B += comp['B_pole']; Barch = np.zeros((N, N))
    for (lo, hi) in ((0.0, L), (L, 2 * L)):
        for sk, wk in zip((hi - lo) / 2 * sg + (lo + hi) / 2, (hi - lo) / 2 * swg):
            ker = -0.5 * (k1(sk) + k2(sk)) + (0.5 * sech(sk) if FORM == 'DH' else 0.0)
            Barch += wk * ker * conv_plus(sk)
    B += Barch; comp['B_arch'] = Barch
    # ---- Toeplitz T ----
    T = np.zeros((N, N))
    if FORM == 'ZETA': T += np.outer(ap, am) + np.outer(am, ap)
    E0 = conv_minus(0.0)   # = identity (orthonormal)
    c0 = -0.5772156649015329 - math.log(math.pi) - math.log(1 - math.exp(-2 * L))
    T += c0 * E0 + (math.log(5) * E0 if FORM == 'DH' else 0)
    for sk, wk in zip(L / 2 * sg + L / 2, L / 2 * swg):
        E = conv_minus(sk)
        T += wk * ((E0 - E) * k1(sk) - E * k2(sk) + (E * sech(sk) if FORM == 'DH' else 0))
    # ---- primes ----
    NN = int(math.exp(2 * L)) + 1
    if FORM == 'ZETA': wts = {n: mang(n) for n in range(2, NN + 1)}
    else: wts = dh_weights(NN)
    Tp = np.zeros((N, N)); Bp = np.zeros((N, N))
    for n in range(2, NN + 1):
        ln = math.log(n)
        if wts.get(n, 0) == 0: continue
        c = wts[n] / math.sqrt(n)
        if ln < 2 * L: Bp += c * conv_plus(ln)
        if ln < L: Tp += 2 * c * conv_minus(ln)
    T -= Tp; B -= Bp; comp['T_prime'] = Tp; comp['B_prime'] = Bp; comp['ap'] = ap; comp['am'] = am
    build.comp = comp
    return 0.5 * (T + T.T), 0.5 * (B + B.T)

if __name__ == '__main__':
    FORM = sys.argv[1]; L = float(eval(sys.argv[2])); N = int(sys.argv[3])
    T, B = build(FORM, L, N)
    eB = np.linalg.eigvalsh(B); eT = np.linalg.eigvalsh(T)
    ePlus = np.linalg.eigvalsh(T + B); eMinus = np.linalg.eigvalsh(T - B)
    out = {"kind": "float", "form": FORM, "L": sys.argv[2], "N": N,
           "B_min": eB[0], "B_max": eB[-1], "B_neg_count": int((eB < -1e-9).sum()), "B_pos_count": int((eB > 1e-9).sum()),
           "T_min": eT[0], "even_min(T+B)": ePlus[0], "odd_min(T-B)": eMinus[0]}
    if eT[0] > 0:
        Ti = np.linalg.cholesky(T); Li = np.linalg.inv(Ti)
        C = Li @ B @ Li.T; eC = np.linalg.eigvalsh(0.5 * (C + C.T))
        out.update({"C_min": eC[0], "C_max": eC[-1], "1-|C|": 1 - max(abs(eC[0]), abs(eC[-1]))})
    print(json.dumps({k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in out.items()}))
