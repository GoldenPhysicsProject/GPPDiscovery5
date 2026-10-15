"""S17 normalization check (NUMERICAL): the assembled Codex/Grok quadratic Q_eta(g) vs the zero side
2*sum_{0<gamma<T} |ghat(gamma)|^2 (real primitive eta: zeros symmetric; ghat(t)=int g e^{itx} dx) using zeros of
L(s,eta) found as sign changes of the real function Lambda(1/2+it,eta) (root number +1 for these real characters).
Tests: conductor sign (+log f), parity gamma constant c_eps, pole only for eta=1, prime sign Re eta(n)."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s17_build import *
import mpmath as mp, numpy as np
mp.mp.dps = 20
TMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 200.0
chars = {"chi3 (mod 3, odd)": (3, 1, [0, 1, -1]), "chi4 (mod 4, odd)": (4, 1, [0, 1, 0, -1]),
         "chi5 quadratic (mod 5, even)": (5, 0, [0, 1, -1, -1, 1])}
def zeros(q, a, chi, T, dt=0.02):
    f = lambda t: (mp.re(mp.power(q/mp.pi, (mp.mpf(1)/2+1j*t+a)/2)*mp.gamma((mp.mpf(1)/2+1j*t+a)/2)*mp.dirichlet(mp.mpf(1)/2+1j*t, chi)))
    ts = np.arange(0.01, T, dt); vals = [f(t) for t in ts]; zs = []
    for i in range(len(ts)-1):
        if vals[i]*vals[i+1] < 0: zs.append(float(mp.findroot(f, (ts[i], ts[i+1]), solver='illinois')))
    return zs
cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"zeros_T{int(TMAX)}.json")
Z = json.load(open(cache)) if os.path.exists(cache) else {}
for nm, (q, a, chi) in chars.items():
    if nm not in Z: Z[nm] = zeros(q, a, chi, TMAX); json.dump(Z, open(cache, "w"))
if "zeta" not in Z:
    Z["zeta"] = []; k = 1
    while True:
        g = float(mp.zetazero(k).imag)
        if g > TMAX: break
        Z["zeta"].append(g); k += 1
    json.dump(Z, open(cache, "w"))
for nm, (q, a, chi) in chars.items():
    T = TMAX; Nexp = T/(2*np.pi)*np.log(q*T/(2*np.pi*np.e)) - (1 if a else -1)/8.0 + 0.0
    print(f"{nm}: {len(Z[nm])} zeros in (0,{T}); smooth count ~ {Nexp:.1f}; first {Z[nm][:3]}")
M = 8
xq, wq = np.polynomial.legendre.leggauss(400)
for hs in ["0.3", "0.5"]:
    B = build(M, hs); F = forms(B, hs); h = float(hs)
    x = xq*h; w = wq*h; u = x/h
    basis = np.array([np.polynomial.legendre.legval(u, [0]*j+[1])*(1-u**2)**4 for j in range(M)])
    for tag, c in (("even g", np.array([1, 0, .3, 0, -.2, 0, 0, 0.])), ("odd g", np.array([0, 1, 0, .4, 0, -.1, 0, 0.])),
                   ("mixed g", np.array([1, .5, .3, -.2, 0, .1, 0, 0.]))):
        g = c @ basis
        for nm in ["zeta"]+list(chars):
            Qm = F[nm]; form = sum(c[i]*c[j]*float(Qm[i, j].mid().str(20, radius=False)) for i in range(M) for j in range(M))
            gam = np.array(Z[nm]); gh = (g*w) @ np.exp(1j*np.outer(x, gam)); zs = 2*np.sum(np.abs(gh)**2)
            print(f"h={hs} {tag:7s} {nm:30s} form={form:+.10e} zero-side(T<{TMAX:.0f})={zs:+.10e} rel.diff={abs(form-zs)/abs(zs):.2e}")
