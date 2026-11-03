"""C.18 shared-entries engine (Arb, python-flint): rigorous enclosures of the band matrices and tail integrals used by
  (i) the Temple-Lehmann route (ARP C.18, Lemma P.17.A applied to the comparison form Q'_T), and
 (ii) Schur-type routes (Thm C.15.A / C.16.A; S&M's Thm S: N_N and n_perp are the same integrals),
for the full Weil form Q_L on real parity-sigma xi in H^1_0(-L,L), general L and prime-power set S (no zero-side input).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN.

Objects (f_k = e_k/n_k H-orthonormal Dirichlet modes, n_k^2 = w_k^2 + 1/4; odd w_k = (k+1)pi/L, even w_k = (k+1/2)pi/L):
  ehat_k(t) = sqrt(L)[sinc((w_k - t)L) + sg sinc((w_k + t)L)], sg = +1 even / -1 odd  (|xi^|; entire)
  NOm_jk  = (1/pi) int_0^T (1/4+t^2)(Omega(T) - Omega(t)) ehat_j ehat_k dt/(n_j n_k)     [= S&M N_N; m_neg Gram]
  Bn_jk   = (1/pi) int_0^T (1/4+t^2)(1 - cos(t log n))    ehat_j ehat_k dt/(n_j n_k)     [band part of the n-line]
  =>  A = q_T I - NOm + sum_n 2c_n Bn   (Gram of the comparison form Q'_T on V),   M = NOm + sum_n 2c_n Bn  (m_abs Gram),
      q_T = Omega(T) - 2 sum_n c_n,  c_n = Lambda(n)/sqrt(n).
  eps^g_K = (1/pi) int_0^T (1/4+t^2) g(t) K_{W_K}(t) dt,  g in {m_abs, m_neg},  K_{W_K} = K_sigma - sum_{k<K} ehat_k^2/n_k^2
  (K_sigma closed form: Lemma C.16.K odd, S&M Lemma A even).
Method (PROVED error bound): [0,T] is cut into panels of width h; on each panel an ng-point Gauss-Legendre rule
(nodes/weights from arb.legendre_p_root, i.e. enclosed). Every integrand f is holomorphic on the closed rectangle
R_p = [c_p - a, c_p + a] x [-b, b] (a = sqrt(b^2 + h^2/4), b = 0.4 < 1/2, the nearest poles of Omega and K_sigma being at
+-i/2), which contains the Bernstein ellipse E_rho, rho = (b + a)/(h/2). If |f| <= M_p on R_p then
|int_p f - G_p(f)| <= (h/2) * 8 M_p rho^(-2 ng)/(rho - 1)   (Chebyshev coefficients |a_k| <= 2 M rho^-k, ATAP Thm 8.1;
Gauss is exact to degree 2ng-1 and has positive weights summing to 2). M_p is bounded by ball evaluation of |f| on boxes
covering the boundary of R_p (maximum modulus), using |ehat_k(z)| <= 2 sqrt(L) e^{L|Im z|} (sinc(u) = (1/2)int_{-1}^1 e^{ius}ds).
The quadrature error radius is added to every entry. All arithmetic is arb/acb ball arithmetic (prec 128 by default).
"""
from flint import arb, acb, arb_mat, ctx
ctx.prec = 128
PI = arb.pi(); LOGPI = PI.log(); Q4 = arb(1)/4

def Omega(z):
    """meromorphic Omega (= Re psi(1/4+it/2) - log pi on the real line); z acb. Poles at +-i(2m+1/2)."""
    h = acb(0, 1)*z/2
    return ((Q4 + h).digamma() + (Q4 - h).digamma())/2 - LOGPI

class Engine:
    def __init__(s, twoL, S, parity, T, K, h=1, ng=40, b='0.4', verbose=True):
        """twoL: arb; S: list of (n, Lambda(n)) with Lambda(n) arb; T: int multiple of h; K: number of modes (rows)."""
        s.twoL = arb(twoL); s.L = s.twoL/2; s.par = parity; s.sg = 1 if parity == 'even' else -1
        s.S = [(n, arb(n).log(), arb(lam)/arb(n).sqrt()) for n, lam in S]  # (n, log n, c_n)
        s.T = int(T); s.h = h; s.ng = ng; s.b = arb(b); s.K = K; s.verbose = verbose
        s.OMT = Omega(acb(s.T)).real; s.QT = s.OMT - 2*sum((c for _, _, c in s.S), arb(0))
        s.w = [((arb(k) + arb(1)/2) if parity == 'even' else arb(k + 1))*PI/s.L for k in range(K)]
        s.n = [(wk*wk + Q4).sqrt() for wk in s.w]
        hh = arb(h)/2; s.a = (s.b**2 + hh**2).sqrt(); s.rho = (s.b + s.a)/hh
        s.qerr_fac = hh*8*s.rho**(-2*ng)/(s.rho - 1)      # times M_p
        s.Ebound = 4*s.L*(2*s.L*s.b).exp()                 # sup |ehat_j ehat_k| on R_p
        roots = [arb.legendre_p_root(ng, i, weight=True) for i in range(ng)]
        s.t = []; s.wt = []
        for p in range(s.T//h):
            c = arb(p*h) + hh
            for x, wx in roots: s.t.append(c + hh*x); s.wt.append(hh*wx)
        s.P = len(s.t)
        s._rows = {}
        s._panel_sup()
    # ---- functions on the real line (arb) and in the complex plane (acb) ----
    def ehat(s, k, t):
        wk = s.w[k]; return s.L.sqrt()*(((wk - t)*s.L).sinc() + s.sg*((wk + t)*s.L).sinc())
    def Kfull(s, t):
        L = s.L; c, sn = (t*L).cos(), (t*L).sin()
        if s.par == 'even': num = L + L*(2*t*L).sinc() - 2*c*((L/2).tanh()/2*c + t*sn)/(t*t + Q4)
        else: num = L - L*(2*t*L).sinc() - 2*sn*(1/(L/2).tanh()/2*sn - t*c)/(t*t + Q4)
        return num/(Q4 + t*t)
    def gOm(s, z): return (Q4 + z*z)*(s.OMT - Omega(z))
    def gn(s, z, lg): return (Q4 + z*z)*(1 - (lg*z).cos())
    def _sup_rect(s, f, c, nb=24):
        """upper bound of |f| on rectangle [c-a,c+a]x[-b,b] (f holomorphic there): boxes covering its boundary."""
        a, b = s.a, s.b; best = arb(0); pts = []
        for i in range(nb):   # horizontal edges
            x = c - a + (2*a)*(arb(2*i + 1)/(2*nb)); r = a/nb
            for y in (b, -b): pts.append((x, y, r, arb(0)))
        nv = max(4, nb//3)
        for i in range(nv):   # vertical edges
            y = -b + (2*b)*(arb(2*i + 1)/(2*nv)); r = b/nv
            for x in (c - a, c + a): pts.append((x, y, arb(0), r))
        for x, y, rx, ry in pts:
            z = acb(arb(x.mid(), rx.upper() if rx != 0 else 0), arb(y.mid(), ry.upper() if ry != 0 else 0))
            v = abs(f(z))
            if not v.is_finite(): raise ValueError('non-finite sup bound')
            u = arb(v.upper()); best = u if u > best else best
        return best
    def _panel_sup(s):
        hh = arb(s.h)/2; s.supOm = []; s.supB = {n: [] for n, _, _ in s.S}; s.supK = []; s.sup4 = []
        for p in range(s.T//s.h):
            c = arb(p*s.h) + hh
            s.supOm.append(s._sup_rect(s.gOm, c))
            for n, lg, _ in s.S: s.supB[n].append(s._sup_rect(lambda z: s.gn(z, lg), c))
            s.supK.append(s._sup_rect(s.Kfull, c))
    # ---- node data ----
    def row(s, k):
        if k not in s._rows: s._rows[k] = [s.ehat(k, t)/s.n[k] for t in s.t]
        return s._rows[k]
    def weights(s, which):
        if which == 'Om': return [w*(Q4 + t*t)*(s.OMT - Omega(acb(t)).real)/PI for t, w in zip(s.t, s.wt)]
        lg = [lg for n, lg, _ in s.S if n == which][0]
        return [w*(Q4 + t*t)*(1 - (lg*t).cos())/PI for t, w in zip(s.t, s.wt)]
    def quad_err(s, which):
        sup = s.supOm if which == 'Om' else s.supB[which]
        return sum((m*s.qerr_fac*s.Ebound/PI for m in sup), arb(0))   # divide by n_j n_k per entry
    def gram(s, which, rows, cols):
        """arb matrix [(1/pi) int_0^T (1/4+t^2) g ehat_j ehat_k/(n_j n_k)] for j in rows, k in cols, g = which-component."""
        wv = s.weights(which); e = s.quad_err(which)
        Fr = arb_mat([[x*w for x, w in zip(s.row(j), wv)] for j in rows])
        Fc = arb_mat([s.row(k) for k in cols]).transpose()
        G = Fr*Fc
        out = [[G[i, l] + arb(0, (e/(s.n[j]*s.n[k])).upper()) for l, k in enumerate(cols)] for i, j in enumerate(rows)]
        return out
    def KW(s, Kp):
        """K_{W_Kp} at the nodes."""
        out = []
        for i, t in enumerate(s.t):
            out.append(s.Kfull(t) - sum((s.row(k)[i]**2 for k in range(Kp)), arb(0)))
        return out
    def eps(s, Kp):
        """(eps^abs_Kp, eps^neg_Kp), enclosed (Gauss sum + quadrature error)."""
        kw = s.KW(Kp); wO = s.weights('Om'); wB = {n: s.weights(n) for n, _, _ in s.S}
        neg = sum((wO[i]*kw[i] for i in range(s.P)), arb(0))
        ab = neg + sum((2*c*sum((wB[n][i]*kw[i] for i in range(s.P)), arb(0)) for n, _, c in s.S), arb(0))
        Ssum = sum((1/(nk*nk) for nk in s.n[:Kp]), arb(0))
        errN = arb(0); errA = arb(0)
        for p in range(s.T//s.h):
            kb = s.supK[p] + s.Ebound*Ssum
            errN += s.supOm[p]*kb*s.qerr_fac/PI
            errA += (s.supOm[p] + sum((2*c*s.supB[n][p] for n, _, c in s.S), arb(0)))*kb*s.qerr_fac/PI
        return ab + arb(0, errA.upper()), neg + arb(0, errN.upper())
    def assemble(s, NOm, B):
        """A (comparison-form Gram) and M (m_abs Gram) from the components; B = {n: matrix}."""
        R = len(NOm); C = len(NOm[0]); A = [[None]*C for _ in range(R)]; M = [[None]*C for _ in range(R)]
        for i in range(R):
            for j in range(C):
                line = sum((2*c*B[n][i][j] for n, _, c in s.S), arb(0))
                M[i][j] = NOm[i][j] + line
                A[i][j] = -NOm[i][j] + line + (s.QT if i == j else 0)
        return A, M

def ldlt_inertia(Mx):
    """ball LDL^T without pivoting; returns (ok, #negative pivots, pivots). ok iff every pivot ball excludes 0
    (then, by Sylvester, #negative pivots = #negative eigenvalues of every symmetric matrix in the box)."""
    n = len(Mx); Lm = [[arb(0)]*n for _ in range(n)]; D = [arb(0)]*n
    for i in range(n):
        for j in range(i + 1):
            sm = Mx[i][j] - sum((Lm[i][k]*Lm[j][k]*D[k] for k in range(j)), arb(0))
            if i == j:
                if not (sm > 0 or sm < 0): return False, None, D[:i] + [sm]
                D[i] = sm; Lm[i][i] = arb(1)
            else: Lm[i][j] = sm/D[j]
    return True, sum(1 for d in D if d < 0), D
