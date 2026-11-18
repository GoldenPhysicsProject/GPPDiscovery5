"""Round 14E (iii): Lehmer-pair probe of a uniform gap on the quotient (Golden Eye falsifier 3).
g(u) = psi(u)[cos(g1 u) - cos(g2 u)], psi = cos^2(pi u/L) on [-L/2, L/2]; (g1, g2) = a pair of zeta ordinates.
14E2 erratum: the zeros are PROOF INPUTS to this falsifier (numerically verified by mp.zetazero): they build the probe and they make
the bound valid (e^(g_i) = 0 for e in the radical BECAUSE g_i are zeros). W(g) itself comes from integer data (Lambda(n), Gamma_R, pole).
Hypothesis: V (resp. the radical vectors e) supported in the window [-L/2, L/2]; E's half-infinite tail is not covered.
Bound (elementary): for e in the radical (e^(g1) = e^(g2) = 0 under RH), |(g-e)^(g1) - (g-e)^(g2)| <= |g1-g2| * int|u||g-e|
  <= |g1-g2| sqrt(L^3/12) ||g-e||, so dist(g, rad) >= |g^(g1) - g^(g2)| / (|g1-g2| sqrt(L^3/12)), and the quotient Rayleigh
  quotient R_q(g) = W(g)/dist(g,rad)^2 <= W(g) (g1-g2)^2 L^3 / (12 |g^(g1)-g^(g2)|^2).   Same bound for any V inside the radical."""
import sys, os, math, json, numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [os.path.join(HERE, '..', 'round14b_euler'), os.path.join(HERE, '..')]
import models14b as md
from scipy.integrate import quad
Z = md.zeta()
def W_direct(g, du, L):
    n = len(g); N = 1 << int(math.ceil(math.log2(2*n)))
    G = np.fft.rfft(g, N); A = np.fft.irfft(np.abs(G)**2, N)[:n]*du     # A(t_k) = int g(x) g(x+t_k) dx, t_k = k du
    t = np.arange(n)*du; A0 = A[0]
    w = Z['w']; arch = Z['C']*A0 + 2*A0*quad(lambda s: float(w(s)), L, np.inf, limit=200)[0]
    f = np.zeros(n); f[1:] = w(t[1:])*(A0 - A[1:])            # integrand ~ (A0-A)/(2t) -> 0 at t=0 (A0 - A = O(t^2))
    arch += 2*np.trapezoid(f, t)
    u = (np.arange(n) - (n-1)/2)*du
    c = np.trapezoid(g*np.cosh(0.5*u), u); s = np.trapezoid(g*np.sinh(0.5*u), u); pole = 2*(c*c - s*s)
    m = np.log(Z['n']) < L; prime = 0.0
    for nn, cc in zip(Z['n'][m], Z['c'][m]):
        prime += 2*cc/math.sqrt(nn)*np.interp(math.log(nn), t, A)
    return arch + pole - prime, dict(arch=arch, pole=pole, prime=prime, norm2=A0)
def ghat(psihat, t, g1, g2): return 0.5*(psihat(t-g1)+psihat(t+g1)-psihat(t-g2)-psihat(t+g2))
def run(L, g1, g2, near):
    du = min(2e-3, 1/(40*g2)); n = int(round(L/du)) + 1; u = np.linspace(-L/2, L/2, n); du = u[1]-u[0]
    psi = np.cos(math.pi*u/L)**2; g = psi*(np.cos(g1*u) - np.cos(g2*u))
    Wg, parts = W_direct(g, du, L)
    def psihat(t):  # int psi(u) e^{itu} du, closed form for cos^2 window
        t = np.asarray(t, float); a = 2*math.pi/L
        with np.errstate(divide='ignore', invalid='ignore'):
            v = np.sin(t*L/2)*(1/t - t/(t*t - a*a))
        return np.where(np.abs(t) < 1e-9, L/2, np.where(np.abs(np.abs(t)-a) < 1e-9, L/4, v))
    dg = abs(ghat(psihat, g1, g1, g2) - ghat(psihat, g2, g1, g2))
    bound = Wg*(g2-g1)**2*L**3/(12*dg**2)
    zsum = float(sum(2*ghat(psihat, z, g1, g2)**2 for z in near))   # +-gamma: check of W = sum |g^(gamma)|^2 (zeros as check)
    return dict(L=L, g1=g1, g2=g2, gap=g2-g1, W=Wg, norm2=parts['norm2'], W_over_norm2=Wg/parts['norm2'], zero_sum_check=zsum,
                dghat=dg, Rq_upper_bound=bound, **{k: parts[k] for k in ['arch', 'pole', 'prime']})
if __name__ == '__main__':
    mp.mp.dps = 20; out = []
    k0 = 6705
    zs = [float(mp.zetazero(k).imag) for k in range(k0, k0+12)]
    i = int(np.argmin(np.diff(zs))); lp = (zs[i], zs[i+1])
    lo = [float(mp.zetazero(k).imag) for k in range(1, 40)]
    for L in [5.0, 7.0, 8.0]:
        r = run(L, *lp, zs); r['pair'] = f'Lehmer #{k0+i},{k0+i+1}'; out.append(r); print(json.dumps(r), flush=True)
        j = i+3 if i+4 < len(zs) else i-3; r = run(L, zs[j], zs[j+1], zs); r['pair'] = f'typical #{k0+j},{k0+j+1}'; out.append(r); print(json.dumps(r), flush=True)
        r = run(L, lo[0], lo[1], lo); r['pair'] = 'first two zeros'; out.append(r); print(json.dumps(r), flush=True)
    json.dump(dict(zeros_near_7005=zs, rows=out), open(os.path.join(HERE, 'lehmer.json'), 'w'), indent=1)

def optimal(L, zs, i, half=20.0):
    """Best probe in the window band [gbar-half, gbar+half] (Dirichlet sines): maximise |g^(g1)-g^(g2)|^2 / W_Z(g),
    W_Z(g) = sum_{zeros in zs} 2|g^(gamma)|^2 (the explicit-formula zero side, CHECKED against the data side below).
    Then R_q <= W(g*) (g2-g1)^2 L^3 / (12 |g*^(g1)-g*^(g2)|^2), with W(g*) recomputed from integer data."""
    g1, g2 = zs[i], zs[i+1]; gb = 0.5*(g1+g2)
    ks = np.arange(max(1, int((gb-half)*L/math.pi)), int((gb+half)*L/math.pi)+1)
    def ph(t):  # int phi_k e^{itu} du, phi_k = sqrt(2/L) sin(k pi (u+L/2)/L)
        k = (ks*math.pi/L)[:, None]; j = ks[:, None]; t = np.atleast_1d(t)[None, :]
        return np.exp(-1j*t*L/2)*math.sqrt(2/L)*k*(1-(-1.0)**j*np.exp(1j*t*L))/(k**2-t**2)
    Z_ = np.array(zs); P = ph(Z_)                     # (nk, nz)
    d = ph(np.array([g1]))[:, 0] - ph(np.array([g2]))[:, 0]   # 14E2 erratum (Golden Eye AUDIT14 6.4): d was undefined (NameError)
    WZ = 2*(P.conj() @ P.T).real                       # real coefficient vectors: |g^(gamma)|^2 = |sum a_k ph_k|^2, +-gamma
    D = np.outer(d.real, d.real) + np.outer(d.imag, d.imag)   # |d.a|^2 for real a
    ev, U = np.linalg.eigh(WZ); keep = ev > 1e-12*ev.max(); Zm = U[:, keep]/np.sqrt(ev[keep])
    mu, V = np.linalg.eigh(Zm.T @ D @ Zm); a = Zm @ V[:, -1]              # maximiser of |d.a|^2 / WZ(a)
    du = 1/(40*g2); n = int(round(L/du))+1; u = np.linspace(-L/2, L/2, n); du = u[1]-u[0]
    g = (a[:, None]*math.sqrt(2/L)*np.sin(np.outer(ks*math.pi/L, u+L/2))).sum(0)
    Wg, parts = W_direct(g, du, L); dg = abs(d @ a)
    return dict(L=L, g1=g1, g2=g2, gap=g2-g1, nbasis=len(ks), nzeros=len(zs), W_data=Wg, W_zero_side=float(a @ WZ @ a), norm2=parts['norm2'],
                dghat=float(dg), q=float(dg**2/Wg), Rq_upper_bound=float(Wg*(g2-g1)**2*L**3/(12*dg**2)), Rayleigh_W_over_norm2=Wg/parts['norm2'],
                ideal_bound=(g2-g1)**2*L**3/24)
