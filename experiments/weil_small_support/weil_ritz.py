"""Small-support Weil form: Ritz minimum over g supported in [-h,h], for zeta, F_theta, DH. NUMERICAL (float).
Basis g_j(x) = P_j(x/h) (1-(x/h)^2)^4, j < M (same family as Claude Code's cert.py). F = g * g~.
Zero-sum functional (explicit formula, no zeros used), E(z) := int g(x) e^{zx} dx:
 zeta-type, applied to F_sig(t) = F(t) e^{sig t}:
   (1/2pi) int w(r) E(ir+sig)E(-ir-sig) dr + E(1/2+sig)E(-1/2-sig) + E(-1/2+sig)E(1/2-sig)
   - sum_n Lambda(n) n^{-1/2} (n^sig + n^{-sig}) F(log n)              [w(r) = Re psi(1/4+ir/2) - log pi]
 zeta: sig=0.  F_theta: sum over sig = +theta, -theta (zeros rho +/- theta).
 DH: (1/2pi) int wD(r)|E(ir)|^2 dr - sum_n 2 c(n) n^{-1/2} F(log n), wD = Re psi(3/4+ir/2) + log(5/pi),
     c(n) = coefficients of -f'/f for f = DH (real coefficients a(n); no pole).
Output: min generalized eigenvalue lam_min(h) = min Q(g)/||g||^2, separately for even and odd g.
Ritz values are UPPER bounds for the true infimum on the window: a negative value is a (float) witness of failure;
a positive value does NOT prove positivity.
Run: python weil_ritz.py M hmin hmax dh
"""
import sys, os, json, numpy as np
from numpy.polynomial import legendre as L
from scipy.special import digamma
from scipy.linalg import eigh
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from common import write_result

M = int(sys.argv[1]) if len(sys.argv) > 1 else 14
hs = np.arange(float(sys.argv[2]), float(sys.argv[3]) + 1e-9, float(sys.argv[4])) if len(sys.argv) > 4 else np.arange(0.2, 1.21, 0.02)
xq0, wq0 = L.leggauss(1200)
kap = (np.sqrt(10 - 2*np.sqrt(5)) - 2) / (np.sqrt(5) - 1)

def vm(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return np.log(p) if m == 1 else 0.0
    return 0.0
def dh_a(n):
    chi = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
    return ((1 - 1j*kap)/2*chi[n % 5] + (1 + 1j*kap)/2*np.conj(chi[n % 5])).real
c_dh = {1: 0.0}
for n in range(2, 400):
    c_dh[n] = dh_a(n)*np.log(n) - sum(c_dh[d]*dh_a(n//d) for d in range(2, n) if n % d == 0)

def basis(x, h):
    u = x / h
    B = np.array([L.legval(u, [0]*j + [1]) for j in range(M)]) * (1 - u**2)**4
    B[:, np.abs(u) >= 1] = 0
    return B
def form(h, kind, theta=0.0):
    xq, wq = xq0*h, wq0*h
    B = basis(xq, h)                                  # M x Q
    gram = (B*wq) @ B.T
    R = 600.0/h; r = np.linspace(0, R, 30001); dr = r[1]-r[0]; wr = np.full_like(r, dr); wr[0] = wr[-1] = dr/2
    E = lambda z: (B*wq) @ np.exp(np.outer(xq, z))   # M x len(z)
    Q = np.zeros((M, M))
    def bil(a, b, weights=None):
        return (a*weights) @ b.T if weights is not None else np.outer(a, b)
    def Fmat(u):   # F_jk(u) = int g_j(x) g_k(x-u) dx, symmetrized
        Bs = basis(xq - u, h); A = (B*wq) @ Bs.T; return (A + A.T)/2
    if kind == "DH":
        w = digamma(0.75 + 0.5j*r).real + np.log(5/np.pi)
        Er = E(1j*r); Q += 2*np.real((Er*(w*wr)) @ np.conj(Er).T)/(2*np.pi)   # r over [0,R], doubled for (-R,0)
        for n in range(2, int(np.exp(2*h)) + 1):
            if c_dh[n] != 0: Q -= 2*c_dh[n]/np.sqrt(n)*Fmat(np.log(n))
    else:
        w = digamma(0.25 + 0.5j*r).real - np.log(np.pi)
        sigs = [0.0] if theta == 0 else [theta, -theta]
        for sg in sigs:
            Ep, Em = E(1j*r + sg), E(-1j*r - sg)
            A = (Ep*(w*wr)) @ Em.T; Q += 2*np.real(A + A.T)/2/(2*np.pi)
            for a in (0.5, -0.5):
                e1, e2 = E(np.array([a + sg]))[:, 0], E(np.array([-a - sg]))[:, 0]
                P = np.outer(e1, e2).real; Q += (P + P.T)/2
            for n in range(2, int(np.exp(2*h)) + 1):
                lv = vm(n)
                if lv: Q -= lv/np.sqrt(n)*(n**sg + n**(-sg))*Fmat(np.log(n))
    return Q, gram
# note on prime normalization: zeta: -sum Lambda n^-1/2 [F(log n)+F(-log n)] = -2 sum ... Fmat (F even).
# F_theta: per sig, -sum Lambda n^-1/2 [F_sig(log n)+F_sig(-log n)] = -sum Lambda n^-1/2 (n^sig+n^-sig) F(log n); summed over 2 sigs.
def lam(h, kind, theta=0.0):
    Q, G = form(h, kind, theta)
    ev, od = list(range(0, M, 2)), list(range(1, M, 2))
    return [float(eigh(Q[np.ix_(I, I)], G[np.ix_(I, I)], eigvals_only=True)[0]) for I in (ev, od)]
cases = [("zeta", "Z", 0.0), ("F_1/4", "Z", 0.25), ("F_1/20", "Z", 0.05), ("DH", "DH", 0.0)]
sel = os.environ.get("CASES"); cases = [c for c in cases if not sel or c[0] in sel.split(",")]
out = {}
for nm, kd, th in cases:
    rows = []
    for h in hs:
        e, o = lam(h, kd, th); rows.append([round(float(h), 4), e, o])
        print(f"{nm:7s} h={h:.3f} even={e:+.6e} odd={o:+.6e}", flush=True)
    neg = [r_[0] for r_ in rows if min(r_[1], r_[2]) < 0]
    out[nm] = {"rows": rows, "first_negative_h_on_grid": neg[0] if neg else None}
tag = f"M{M}_{hs[0]:.3f}_{hs[-1]:.3f}" + ("_" + sel.replace("/","").replace(",","_") if sel else "")
write_result(os.path.join(os.path.dirname(__file__), f"result_{tag}.json"),
  {"experiment": "weil_small_support", "label": "NUMERICAL (Ritz upper bounds)", "M": M, "results": out},
  kind="float", method="Legendre-bump Ritz, Gauss-Legendre 1200 + trapezoid in r")
