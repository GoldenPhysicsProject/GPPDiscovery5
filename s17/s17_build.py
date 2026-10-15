"""S17 Arb builder: completed primitive-Dirichlet Weil form on C_h = C_c^inf((-h,h);R), Legendre*(1-u^2)^4 basis
(reuses S14 exact rational autocorrelation polynomials G_jk: F_jk(x) = h G_jk(x/h), x in [0,2h]).
Codex normalization (work/codex/2026-10-15_exact_completed_dirichlet_family_weil_K_operator.md):
  Q_eta = 1_{eta=1} P + E_eps + (log f - c_eps) N - Wprime_eta,
  E_eps = int_0^inf nu_eps(t)(F(0)-F(t)) dt = K_eps + T_eps(2h) N   (F=0 beyond 2h),
  nu_eps(t)=2e^{-(1+2eps)t/2}/(1-e^{-2t}),  c_eps = log pi - psi((1+2eps)/4),
  P = int F(t) 2cosh(t/2) dt = 2 M(1/2)M(-1/2),
  Wprime_eta = 2 sum_n Lambda(n) Re eta(n) n^{-1/2} F(log n)   (real g).
F_theta control (S4 weil_ritz.py convention: zeros rho +/- theta, summed over sigma=+/-theta):
  Q_theta = 2[K_0 - Kth + T_0(2h) N - c_0 N] + Pth - 2 sum Lambda(n) n^{-1/2}(n^th+n^-th) F(log n),
  Kth = int_0^{2h} nu_0 F (cosh(theta t)-1) dt,  Pth = 2 int F(t) 2cosh(t/2) cosh(theta t) dt.
All moments int_0^2 v^m phi(v) dv by acb.integral (rigorous enclosures)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s14"))
from s14_exact import G, legendre, p
from flint import fmpq, arb, acb, ctx, arb_mat
ctx.prec = 256
_GC = {}
def Gs_for(M):
    if M not in _GC:
        Gs = {}
        for j in range(M):
            for k in range(j, M):
                if (j+k) % 2 == 0:
                    g = (G(j, k)+G(k, j))/2; Gs[j, k] = Gs[k, j] = g
        _GC[M] = Gs
    return _GC[M]
def vM(n):
    for q in range(2, n+1):
        if n % q == 0:
            m = n
            while m % q == 0: m //= q
            return q if m == 1 else 0
    return 0
def build(M, h, theta=fmpq(1, 4), prec=256, nmax_primes=None):
    """returns dict of arb_mat: N, K0, K1, P, Kth, Pth, and T[n] for prime powers n with log n < 2h."""
    ctx.prec = prec
    h = arb(h) if not isinstance(h, arb) else h
    Gs = Gs_for(M); D = max(g.degree() for g in Gs.values())+1
    hc = acb(h); I = acb(0, 1); th = acb(theta)
    def mom(phi):
        return [acb.integral(lambda z, an, m=m: z**m*phi(z), 0, 2).real for m in range(D+1)]
    xnu0 = lambda z: (hc*z/2).exp()/(I*hc*z).sinc()          # x nu_0(x), x = h z
    xnu1 = lambda z: (-hc*z/2).exp()/(I*hc*z).sinc()         # x nu_1(x)
    mK0, mK1 = mom(xnu0), mom(xnu1)
    mP = mom(lambda z: (hc*z/2).cosh())
    mPth = mom(lambda z: (hc*z/2).cosh()*(th*hc*z).cosh())
    # nu_0(x)(cosh(th x)-1) = [x nu_0(x)] * (th^2 x/2) sinc(i th x/2)^2 ;  x = h z
    mKth = mom(lambda z: xnu0(z)*(th*th*hc*z/2)*((I*th*hc*z/2).sinc())**2)
    ln = {}
    n = 2
    while True:
        L = arb(n).log()
        if not (L < 2*h): break
        if vM(n): ln[n] = L/h
        n += 1
    out = {k: arb_mat(M, M) for k in ("N", "K0", "K1", "P", "Kth", "Pth")}
    T = {n: arb_mat(M, M) for n in ln}
    for (j, k), g in Gs.items():
        co = [arb(c.p)/c.q for c in g.coeffs()]
        q = [-c for c in co[1:]]
        out["N"][j, k] = h*co[0]
        out["K0"][j, k] = h*sum(q[m]*mK0[m] for m in range(len(q)))
        out["K1"][j, k] = h*sum(q[m]*mK1[m] for m in range(len(q)))
        out["P"][j, k] = 4*h*h*sum(co[m]*mP[m] for m in range(len(co)))
        out["Pth"][j, k] = 8*h*h*sum(co[m]*mPth[m] for m in range(len(co)))
        out["Kth"][j, k] = h*h*sum(co[m]*mKth[m] for m in range(len(co)))
        for n, v in ln.items():
            T[n][j, k] = h*sum(co[m]*v**m for m in range(len(co)))
    out["T"] = T
    return out
def Tarb(L, e):
    r = (-arb(L)/2).exp(); at = ((1+r)/(1-r)).log()/2; an = r.atan()
    return 2*(at+an) if e == 0 else 2*(at-an)
def consts():
    C0 = (8*arb.pi()).log()+arb.const_euler()+arb.pi()/2
    return C0, C0-arb.pi()
# characters: (name, conductor, parity, Re eta(n) for n coprime handled by function)
def chi3(n): return [0, 1, -1][n % 3]
def chi4(n): return [0, 1, 0, -1][n % 4]
def chi5q(n): return [0, 1, -1, -1, 1][n % 5]                 # Legendre (n/5), even
def chi5c_re(n): return [0, 1, 0, 0, -1][n % 5]               # Re of quartic chi mod 5 (chi(2)=i), odd
CHARS = [("chi3 (mod 3, odd)", 3, 1, chi3), ("chi4 (mod 4, odd)", 4, 1, chi4),
         ("chi5 quadratic (mod 5, even)", 5, 0, chi5q), ("chi5 quartic (mod 5, odd; = conj on real g)", 5, 1, chi5c_re)]
def forms(B, h, theta=fmpq(1, 4)):
    """assemble arb_mat Q for zeta, F_theta, and the four primitive characters."""
    h = arb(h) if not isinstance(h, arb) else h
    C0, C1 = consts(); T0, T1 = Tarb(2*h, 0), Tarb(2*h, 1)
    N = B["N"]; res = {}
    W = lambda wfun: sum((B["T"][n]*(2*arb(vM(n)).log()*wfun(n)/arb(n).sqrt()) for n in B["T"]), arb_mat(N.nrows(), N.ncols()))
    res["zeta"] = B["P"]+B["K0"]+N*(T0-C0)-W(lambda n: 1)
    tha = arb(theta)
    res["F_theta"] = (B["K0"]-B["Kth"]+N*(T0-C0))*2+B["Pth"]-W(lambda n: (arb(n)**tha+arb(n)**(-tha)))
    for nm, f, e, ch in CHARS:
        Ke, Te, Ce = (B["K0"], T0, C0) if e == 0 else (B["K1"], T1, C1)
        res[nm] = Ke+N*(Te+arb(f).log()-Ce)-W(lambda n, ch=ch: ch(n))
    return res
