"""Grok S97 -- H96a: discrete Krein/de Branges chain for the nested window family.
Nested basis: step functions chi_k on cells [kh,(k+1)h], k=0..n-1 (window [0,T], T=nh). Weil form is translation
invariant, so the Gram is Toeplitz: c_m = W(E_m), E_m = even part of tent_h(z - m h). Levinson-Durbin gives the
Schur/Verblunsky (reflection) coefficients kappa_k and pivots P_k = P_{k-1}(1-kappa_k^2): T_n PD <=> c_0>0, |kappa_k|<1.
In the h->0 limit kappa_k/h -> Krein's continuous reflection function, which fixes the de Branges Hamiltonian H(t).
Distribution (calibrated against s87/hp_forms at N=1, see calib()):
 W(g) = c00 g(0) + int_0^inf k1(z)(g(0)-w(z) g(z)) dz - int_0^inf w(z) k2'(z) g(z) dz - 2 sum_n a_n w(log n)/sqrt(n) g(log n) + poles,
 g even, w = cosh(theta z) (F_theta), k2' = k2 - sech (DH), c00 = -gamma - log pi (+log 5 for DH).
HIGH-PRECISION NUMERICAL (mpmath), not interval. RH OPEN."""
import sys, mpmath as mp
sys.path.insert(0, "../s87"); import hp_forms as H

def setup(dps=60): mp.mp.dps = dps

def spec(FORM, theta=0, NN=10**4, drop=()):
    th = mp.mpf(theta); dh = FORM == "DH"
    w = (lambda z: mp.cosh(th * z)) if FORM == "FT" else (lambda z: mp.mpf(1))
    k1 = lambda z: mp.exp(-z) / mp.sinh(z)
    k2 = lambda z: 2 * mp.exp(-z / 4) * mp.sinh(3 * z / 4) / mp.sinh(z) - (1 / mp.cosh(z / 2) if dh else 0)
    if dh: a = H.dh_weights(NN); a = {n: a[n] for n in a if n >= 2 and a[n] != 0}
    else: a = {n: H.mang(n) for n in range(2, NN + 1) if H.mang(n) != 0}
    for p in drop:
        for n in list(a):
            m = n
            while m % p == 0: m //= p
            if m == 1: del a[n]
    poles = {"ZETA": [(mp.mpf(1) / 2, 2)], "FT": [(mp.mpf(1) / 2 + th, 1), (mp.mpf(1) / 2 - th, 1)]}.get(FORM, [])
    return dict(w=w, k1=k1, k2=k2, a=a, poles=poles, c00=-mp.euler - mp.log(mp.pi) + (mp.log(5) if dh else 0))

def W_tent(S, h, m, pole_scale=1):
    """W applied to g = even part of tent_h(z - m h) (tent: max(0, h - |z|))."""
    h = mp.mpf(h)
    lo, hi = max(0, (m - 1)) * h, (m + 1) * h
    if m == 0:
        g = lambda z: max(h - abs(z), 0)
        g0 = h
    else:
        g = lambda z: (max(h - abs(z - m * h), 0) + max(h - abs(z + m * h), 0)) / 2
        g0 = mp.mpf(0)
    w, k1, k2 = S["w"], S["k1"], S["k2"]
    pts = [lo] + ([m * h] if m > 0 else []) + [hi]
    val = S["c00"] * g0
    if m == 0:
        val += mp.quad(lambda z: k1(z) * (g0 - w(z) * g(z)), pts) + g0 * mp.quad(k1, [hi, mp.inf])
    else:
        val += -mp.quad(lambda z: k1(z) * w(z) * g(z), pts)
    val += -mp.quad(lambda z: w(z) * k2(z) * g(z), pts)
    for n, an in S["a"].items():
        ln = mp.log(n)
        if ln >= hi: break
        if ln > lo: val -= 2 * an * w(ln) / mp.sqrt(n) * g(ln)
    # poles: sum_poles c * (ghat(al) + ghat(-al))/2-type term; ghat(al) = int g(z) e^{al z} dz over R (g even)
    for al, c in S["poles"]:
        gh = 2 * mp.quad(lambda z: g(z) * mp.cosh(al * z), pts)
        val += pole_scale * c * gh
    return val

def toeplitz_c(S, h, n, pole_scale=1):
    return [W_tent(S, h, m, pole_scale) for m in range(n)]

def levinson(c):
    """Durbin on symmetric Toeplitz c[0..n-1]. Returns kappas, pivots P_k (prediction errors), first k where PD fails."""
    n = len(c); P = [c[0]]; kap = []; a = []
    fail = None if c[0] > 0 else 0
    for k in range(1, n):
        acc = c[k] - mp.fsum(a[j] * c[k - 1 - j] for j in range(len(a)))
        kk = acc / P[-1]
        a = [a[j] - kk * a[len(a) - 1 - j] for j in range(len(a))] + [kk]
        kap.append(kk); P.append(P[-1] * (1 - kk ** 2))
        if fail is None and P[-1] <= 0: fail = k
    return kap, P, fail

def calib(L=0.5):
    """Compare n=2 cells of width L (even/odd extension of a constant on [0,L]) with hp_forms N=1 Q_e, Q_o."""
    out = {}
    for F, th in (("ZETA", 0), ("DH", 0), ("FT", 0.1)):
        d = H.build(F, L, 1, th)
        Qe = d["Ae"][0][0] - d["Pe"][0][0] + sum(v[0] ** 2 for v in d["pe"])
        Qo = d["Ao"][0][0] - d["Po"][0][0] - sum(v[0] ** 2 for v in d["po"])
        S = spec(F, th); c = toeplitz_c(S, L, 2)
        e = (c[0] + c[1]) * 2 / L; o = (c[0] - c[1]) * 2 / L   # v=(1,1)/sqrt L -> 2(c0+c1)/L
        out[F] = (mp.nstr(Qe, 15), mp.nstr(e, 15), mp.nstr(Qo, 15), mp.nstr(o, 15))
    return out
if __name__ == "__main__":
    setup(30)
    for k, v in calib().items(): print(k, v)
