"""New Bot (rep-theory seat): GL(1) twist calculus test. NUMERICAL (float64). RH OPEN.
Builds the S14-normalised window Weil form for kernel g(r)*W_zeta(r) where
  g = 1 (zeta), g = 2cosh(theta r) (F_theta = Lam(s-th)Lam(s+th)), g = 2cos(th' r) (G = Lam(s+i th')Lam(s-i th')).
Adapted (copied, not modified in place) from GPPDiscovery5/s82os/os_hankel.py (Grok S82-OS), same conventions.
"""
import math, json, sys
import numpy as np
from numpy.polynomial import legendre as leg
sys.path.insert(0, '/workspace/GPPDiscovery5/s82os')
from os_hankel import mang

def build(L, N, g, poles, M=96, S=160):
    # poles: list of (c_a, a) with 2cosh(r/2) g(r) = sum c_a exp(a r)
    xg, wg = leg.leggauss(M)
    def basis(x):
        x = np.asarray(x, float); t = 2 * x / L - 1
        out = np.zeros((N,) + x.shape); inside = (x >= 0) & (x <= L)
        for k in range(N):
            c = np.zeros(k + 1); c[k] = 1
            out[k] = np.where(inside, math.sqrt((2 * k + 1) / L) * leg.legval(t, c), 0.0)
        return out
    def conv_plus(s):
        a, b = max(0.0, s - L), min(L, s)
        if b <= a: return np.zeros((N, N))
        x = (b - a) / 2 * xg + (a + b) / 2; w = (b - a) / 2 * wg
        return (basis(x) * w) @ basis(s - x).T
    def conv_minus(r):
        a, b = r, L
        if b <= a: return np.zeros((N, N))
        x = (b - a) / 2 * xg + (a + b) / 2; w = (b - a) / 2 * wg
        Cm = (basis(x) * w) @ basis(x - r).T
        return 0.5 * (Cm + Cm.T)
    k1 = lambda z: math.exp(-z) / math.sinh(z)
    k2 = lambda z: 2 * math.exp(-z / 4) * math.sinh(3 * z / 4) / math.sinh(z)
    x = L / 2 * xg + L / 2; w = L / 2 * wg; P = basis(x)
    mom = lambda a: P @ (w * np.exp(a * x))
    B = np.zeros((N, N), complex); T = np.zeros((N, N), complex)
    for c, a in poles:
        B += c * np.outer(mom(a), mom(a)); T += c * np.outer(mom(a), mom(-a))
    B = B.real; T = T.real
    sg, swg = leg.leggauss(S)
    for (lo, hi) in ((0.0, L), (L, 2 * L)):
        for sk, wk in zip((hi - lo) / 2 * sg + (lo + hi) / 2, (hi - lo) / 2 * swg):
            B += wk * (-0.5 * (k1(sk) + k2(sk))) * g(sk) * conv_plus(sk)
    E0 = conv_minus(0.0)
    c0 = -0.5772156649015329 - math.log(math.pi) - math.log(1 - math.exp(-2 * L))
    T += c0 * g(0.0) * E0
    for sk, wk in zip(L / 2 * sg + L / 2, L / 2 * swg):
        E = conv_minus(sk)
        T += wk * ((g(0.0) * E0 - g(sk) * E) * k1(sk) - g(sk) * E * k2(sk))
    NN = int(math.exp(2 * L)) + 1
    for n in range(2, NN + 1):
        lam = mang(n)
        if lam == 0: continue
        ln = math.log(n); c = lam / math.sqrt(n) * g(ln)
        if ln < 2 * L: B -= c * conv_plus(ln)
        if ln < L: T -= 2 * c * conv_minus(ln)
    T = 0.5 * (T + T.T); B = 0.5 * (B + B.T)
    return T + B, T - B, basis

def forms(kind, L, N, th):
    if kind == 'ZETA':
        return build(L, N, lambda r: 1.0, [(1, .5), (1, -.5)])
    if kind == 'F':
        return build(L, N, lambda r: 2 * math.cosh(th * r), [(1, .5 + th), (1, -.5 - th), (1, .5 - th), (1, -.5 + th)])
    if kind == 'G':
        return build(L, N, lambda r: 2 * math.cos(th * r), [(1, .5 + 1j * th), (1, .5 - 1j * th), (1, -.5 + 1j * th), (1, -.5 - 1j * th)])

if __name__ == '__main__':
    N = 24; out = []
    xs = np.linspace(0, 1, 4001)
    for L in [0.3, 0.4, 0.5, 0.6, 0.7]:
        Ze, Zo, basis = forms('ZETA', L, N, 0)
        ez, vz = np.linalg.eigh(Ze); ezo, vzo = np.linalg.eigh(Zo)
        row = {"L": L, "zeta_even_min": ez[0], "zeta_odd_min": ezo[0]}
        for th in [0.05, 0.1, 0.25]:
            Fe, Fo, _ = forms('F', L, N, th); Ge, Go, _ = forms('G', L, N, th)
            ef, vf = np.linalg.eigh(Fe); eg = np.linalg.eigvalsh(Ge); ego = np.linalg.eigvalsh(Go); efo = np.linalg.eigvalsh(Fo)
            u = xs * L; Bu = basis(u)
            fstar = vf[:, 0] @ Bu; zstar = vz[:, 0] @ Bu
            wv = fstar * np.cosh(th * u)
            ov = abs(np.trapezoid(wv * zstar, u)) / math.sqrt(np.trapezoid(wv * wv, u) * np.trapezoid(zstar * zstar, u))
            nodes = int(np.sum(np.diff(np.sign(fstar[np.abs(fstar) > 1e-12])) != 0))
            row[f"th{th}"] = {"F_even_min": ef[0], "F_odd_min": efo[0], "G_even_min": eg[0], "G_odd_min": ego[0],
                              "overlap(F_vec*cosh, zeta_even_vec)": ov, "F_vec_nodes_half": nodes}
        out.append(row); print(json.dumps(row, default=float))
    json.dump(out, open('twist_results.json', 'w'), indent=1, default=float)
