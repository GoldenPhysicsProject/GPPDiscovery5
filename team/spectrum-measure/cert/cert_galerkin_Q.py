"""CERT (Arb ball arithmetic, python-flint): sector coercivity of the p = 2 Weil form Q_L = R - 2c2 B on H^1_0(-L, L),
2L = log 2 + ell, by a capped-symbol Galerkin/Schur argument with rigorous tail (CERT_p2_sectors.md).
Spectrum & Measure seat, GPP RH team. Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN.

Steps (all quantities arb/acb balls; integrals by acb.integral on holomorphic integrands; non-finite => abort):
  Om(t) = (psi(1/4+it/2)+psi(1/4-it/2))/2 - log pi  (= Re psi(1/4+it/2) - log pi on the real line), Om* = Om(T2).
  R >= Om* ||.||_H^2 - N,  N(xi) = (1/pi) int_0^T2 (Om* - Om)(1/4+t^2)|xi^|^2   (Om increasing, Prime (Om2)).
  V_N = span of the first N H-orthonormal Dirichlet modes f_m of the sector; xi = u + w, w H-orthogonal to V_N.
  N(xi) <= (1+r) N(u) + (1+1/r) n_perp ||w||^2, n_perp = int_0^T2 (Om*-Om) dens_N, dens_N = dens - (1/pi)(1/4+t^2) sum_m f_m^2 (Lemma A).
  P := Om*||.||^2 - 2c2 B;  |B(u,w)| <= (1/2)||u||_{H(E)} ||w||_{H(E)},  |B(w)| <= (1/2)||w||_{H(E)}^2 <= (1/2)||w||^2.
  => Q(xi) >= u^T M u + eta*beta0 ||w||^2, M = Om* I - (1+r) N_N - 2c2 B_N - (c2^2/beta) G_E,
     beta0 = Om* - c2 - (1+1/r) n_perp, beta = (1-eta) beta0.  M >= lam0 certified by ball Cholesky of M - lam0 I.
  Margin: Q(xi) >= min(lam0, eta*beta0) ||xi||_H^2 for all xi in the sector of H^1_0(-L,L); monotone in ell (form L-independent).
"""
import sys, time
from flint import arb, acb, ctx, arb_mat

def chk(x, what):
    ok = x.is_finite() if isinstance(x, arb) else (x.real.is_finite() and x.imag.is_finite())
    if not ok: raise SystemExit("NON-FINITE BALL: %s = %s" % (what, x))
    return x

def run(par, ell, N, T2, r, eta, lam0):
    pi = arb.pi(); a = arb(2).log(); c2 = a/arb(2).sqrt(); q = arb(1)/4
    ell = arb(ell); L = (a + ell)/2; T2i = int(T2); T2 = arb(T2); r = arb(r); eta = arb(eta); lam0 = arb(lam0)
    I = acb(0, 1); lp = acb(pi.log())
    def Om(z): return ((acb(q) + I*z/2).digamma() + (acb(q) - I*z/2).digamma())/2 - lp
    Os = chk(Om(acb(T2)).real, "Om*")
    sg = 1 if par == 'e' else -1
    ks = [((2*m+1) if par == 'e' else 2*(m+1))*pi/(2*L) for m in range(N)]
    nr = [1/(L*(k*k + q)).sqrt() for k in ks]
    def fh(z, m): k = ks[m]; return (((z-k)*L).sinc() + sg*((z+k)*L).sinc())*L*nr[m]
    # N_N
    NN = [[None]*N for _ in range(N)]
    for i in range(N):
        for j in range(i, N):
            f = lambda z, an, i=i, j=j: (Os - Om(z))*(acb(q) + z*z)*fh(z, i)*fh(z, j)
            v = sum((acb.integral(f, s, s+10) for s in range(0, T2i, 10)), acb(0))/pi
            v = chk(v, "N[%d,%d]" % (i, j)).real
            NN[i][j] = NN[j][i] = v
    # n_perp
    if par == 'e':
        dens = lambda z: (acb(L) + L*(2*z*L).sinc() - 2*(z*L).cos()*((L/2).tanh()/2*(z*L).cos() + z*(z*L).sin())/(z*z + q))/pi
    else:
        dens = lambda z: (acb(L) - L*(2*z*L).sinc() - 2*(z*L).sin()*((L/2).cosh()/(L/2).sinh()/2*(z*L).sin() - z*(z*L).cos())/(z*z + q))/pi
    def fp_(z, an):
        hs = acb(0)
        for m in range(N): hs += fh(z, m)**2
        return (Os - Om(z))*(dens(z) - (acb(q) + z*z)*hs/pi)
    nperp = chk(sum((acb.integral(fp_, s, s+5) for s in range(0, T2i, 5)), acb(0)), "nperp").real
    # B_N and G_E on S- = [-L, -L+ell]
    if par == 'e':
        f0 = lambda z, m: (ks[m]*z).cos()*nr[m]; f1 = lambda z, m: -ks[m]*(ks[m]*z).sin()*nr[m]
    else:
        f0 = lambda z, m: (ks[m]*z).sin()*nr[m]; f1 = lambda z, m: ks[m]*(ks[m]*z).cos()*nr[m]
    lo = -L; hi = -L + ell; aa = acb(a)
    BN = [[None]*N for _ in range(N)]; GE = [[None]*N for _ in range(N)]
    for i in range(N):
        for j in range(i, N):
            if ell == 0:
                b = arb(0); g = arb(0)
            else:
                fb = lambda z, an, i=i, j=j: (f1(z, i)*f1(z+aa, j) + f1(z, j)*f1(z+aa, i) + q*(f0(z, i)*f0(z+aa, j) + f0(z, j)*f0(z+aa, i)))/2
                fg = lambda z, an, i=i, j=j: 2*(f1(z, i)*f1(z, j) + q*f0(z, i)*f0(z, j))
                b = chk(acb.integral(fb, lo, hi), "B").real; g = chk(acb.integral(fg, lo, hi), "G").real
            BN[i][j] = BN[j][i] = b; GE[i][j] = GE[j][i] = g
    beta0 = Os - c2 - (1 + 1/r)*nperp; beta = (1 - eta)*beta0
    M = [[(Os if i == j else arb(0)) - (1+r)*NN[i][j] - 2*c2*BN[i][j] - c2*c2/beta*GE[i][j] - (lam0 if i == j else arb(0))
          for j in range(N)] for i in range(N)]
    # ball Cholesky (LDL^T): all pivots must be certified > 0
    A = [row[:] for row in M]; minpiv = None
    for k in range(N):
        p = A[k][k]
        if not p > 0: return dict(ok=False, k=k, piv=p, Os=Os, nperp=nperp, beta0=beta0)
        minpiv = p if minpiv is None or p.upper() < minpiv.upper() else minpiv
        for i in range(k+1, N):
            lik = A[i][k]/p
            for j in range(k+1, i+1):
                A[i][j] = A[i][j] - lik*A[j][k]; A[j][i] = A[i][j]
    return dict(ok=True, Os=Os, nperp=nperp, beta0=beta0, wmargin=eta*beta0, lam0=lam0, minpiv=minpiv)

if __name__ == "__main__":
    ctx.prec = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    print("prec bits:", ctx.prec)
    jobs = [("e", "0", 40, 100, "0.003", "0.5", "0.05"),
            ("e", "0.06", 40, 100, "0.003", "0.1", "0.009"),
            ("o", "0.15", 40, 100, "0.003", "0.1", "0.06")]
    if len(sys.argv) > 2:   # single job from the command line: par ell N T2 r eta lam0
        p = sys.argv[2:]; jobs = [(p[0], p[1], int(p[2]), int(p[3]), p[4], p[5], p[6])]
    for jb in jobs:
        t0 = time.time(); res = run(*jb)
        print("sector=%s ell=%s N=%d T2=%d r=%s eta=%s lam0=%s ->" % jb, res, "(%.0fs)" % (time.time()-t0)); sys.stdout.flush()
