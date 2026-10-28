"""C.15 Arb callbacks: acb.integral-ready integrands for the even p=2 certificate (Thm C.15.A) at ell=0.1, T=60, N=64.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN.

Every callback has the python-flint signature f(x, analytic) -> acb, for acb.integral(f, a, b).

Quantities (all on [0, T]; see ARP section C.15):
  A_jk  = q_T delta_jk + (1/pi) int_0^T (1/4+t^2)(q(t)-q_T) ehat_j ehat_k dt / (n_j n_k)        -> mu_A (interval Cholesky of A - mu I)
  M_jk  =               (1/pi) int_0^T (1/4+t^2) m_abs(t)   ehat_j ehat_k dt / (n_j n_k)        -> s_V (Cholesky of s I - M, or Frobenius)
  eps^g = (1/pi) int_0^T (1/4+t^2) g(t) K_W(t) dt,  g in {m_abs, m_neg},  K_W = K_e - sum_{k<N} ehat_k^2/n_k^2
  q(t) = Omega(t) - 2 c2 cos(a t), q_T = Omega(T) - 2 c2, m_abs = Omega(T)-Omega(t)+2c2(1-cos at), m_neg = Omega(T)-Omega(t),
  Omega(t) = (psi(1/4+it/2) + psi(1/4-it/2))/2 - log pi,  ehat_k(t) = sqrt(L)[sinc((w_k-t)L) + sinc((w_k+t)L)],  n_k^2 = w_k^2+1/4.

DOMAIN / ANALYTICITY CHOICES (documented per Golden Eye's Arb requirements):
 (1) Omega: we NEVER use Re psi (non-holomorphic). We use the symmetrised form (psi(1/4+it/2)+psi(1/4-it/2))/2, which equals
     Re psi on the real line and is MEROMORPHIC in t, with poles only at t = +-i(2m+1/2), m>=0, i.e. |Im t| >= 1/2.
     acb.integral's Gauss-Legendre ellipses around real subintervals hit those poles only as non-finite balls, which acb.integral
     handles by bisection; no branch cut exists, so the analytic flag may be ignored (flint docs: meromorphic case).
 (2) ehat_k: written with acb.sinc (entire); the removable singularity at t = w_k of 2 w cos(tL)/(w^2-t^2) never appears.
 (3) K_e (S&M Lemma A): sin(2tL)/(2t) is written as L*sinc(2tL) (entire); the factor 1/(t^2+1/4) has poles at +-i/2 only
     (meromorphic). The closed form K_e - sum_{k<N} cancels to relative size ~1e-5: run at prec >= 128 bits.
 (4) Majorants: the RECOMMENDED m_abs, m_neg are analytic (no abs, no max), so every recommended integrand is meromorphic.
     The sharp choices |a_T| and (a_T)_- are NOT holomorphic; their callbacks return a non-finite ball when analytic=True
     (forcing acb.integral to refine to point evaluations, slow and only O(h) accurate). To use them rigorously, split [0,T] at
     Arb-isolated zeros of a_T and integrate +-a_T (analytic) on each piece. Not needed: the analytic majorants already give
     float lambda_cert = 0.01385 at ell = 0.1.
 (5) |t| >= T: NOTHING is integrated there. Step 1 of Thm C.15.A only uses q(t) - q_T >= 0 on |t| >= T, which is PROVED
     analytically (Omega increasing: d/dy Re psi(x+iy) = sum 2y/((n+x)^2+y^2)^2 > 0; 1 - cos >= 0). The majorants live on [0,T].
 (6) Splitting: we integrate [0,T] in panels of width 5 (PANEL) purely for speed; every panel is a real segment, no singularity
     lies on it. Tail for the mode-sum alternative to (3): sum_{k>=K} f_k^2 <= (L/pi)*4/(3 L (w_K - pi/L - T)^3) (C.15.K).
 Labels: the callbacks and splitting are PROVED-correct as formulas; any number printed by the self-test is NUMERICAL until
 S&M's full run encloses every input and certifies mu_A / s_V by interval Cholesky (REMARK: Arb-run hypotheses of C.15.D).
"""
from flint import acb, arb, ctx
ctx.prec = 128
ELL = arb('0.1'); A2 = arb(2).log(); C2 = A2/arb(2).sqrt(); L = (A2 + ELL)/2; T = arb(60); N = 64; PANEL = 5
PI = arb.pi(); LOGPI = PI.log(); Q4 = arb(1)/4

def Omega(t):
    """meromorphic Omega(t); t acb."""
    h = acb(0, 1)*t/2
    return ((Q4 + h).digamma() + (Q4 - h).digamma())/2 - LOGPI
OMT = Omega(acb(T)).real
QT = OMT - 2*C2
def w(k): return (arb(k) + arb(1)/2)*PI/L
def nk(k): return (w(k)**2 + Q4).sqrt()
def ehat(k, t):
    wk = w(k); return L.sqrt()*(((wk - t)*L).sinc() + ((wk + t)*L).sinc())
def q(t): return Omega(t) - 2*C2*(A2*t).cos()
def m_abs(t): return OMT - Omega(t) + 2*C2*(1 - (A2*t).cos())
def m_neg(t): return OMT - Omega(t)

def A_integrand(j, k):
    """integrand of (A_jk - q_T delta_jk); meromorphic (choice 1,2)."""
    s = 1/(PI*nk(j)*nk(k))
    return lambda t, analytic: s*(Q4 + t*t)*(q(t) - QT)*ehat(j, t)*ehat(k, t)
def M_integrand(j, k):
    """integrand of M_jk (m_abs-weighted band Gram matrix, s_V = lambda_max(M)); meromorphic (choice 4)."""
    s = 1/(PI*nk(j)*nk(k))
    return lambda t, analytic: s*(Q4 + t*t)*m_abs(t)*ehat(j, t)*ehat(k, t)
def Ke(t):
    """S&M Lemma A even kernel sum_{k>=0} ehat_k^2/n_k^2; meromorphic (choice 3)."""
    c, s = (t*L).cos(), (t*L).sin()
    num = L + L*(2*t*L).sinc() - 2*c*((L/2).tanh()/2*c + t*s)/(t*t + Q4)
    return num/(Q4 + t*t)
def KW(t, Nm=N):
    return Ke(t) - sum(ehat(k, t)**2/(w(k)**2 + Q4) for k in range(Nm))
def eps_integrand(g, Nm=N):
    """integrand of eps^g_N, g in {m_abs, m_neg}; meromorphic."""
    return lambda t, analytic: (Q4 + t*t)*g(t)*KW(t, Nm)/PI
def sharp_abs_integrand(j, k):
    """NOT holomorphic (|a_T|): returns non-finite when analytic=True (choice 4). Use only with zero-splitting."""
    s = 1/(PI*nk(j)*nk(k))
    def f(t, analytic):
        if analytic: return acb('nan')
        return s*(Q4 + t*t)*acb(abs((q(t) - QT).real))*ehat(j, t)*ehat(k, t)
    return f
def integrate(f, a=0, b=None, panel=PANEL, **kw):
    b = T if b is None else arb(b); a = arb(a); tot = acb(0); x = a
    while x < b:
        y = x + panel
        if y > b: y = b
        tot += acb.integral(f, x, y, **kw); x = y
    return tot
def A_entry(j, k): return integrate(A_integrand(j, k)).real + (QT if j == k else 0)
def M_entry(j, k): return integrate(M_integrand(j, k)).real
def eps_entry(g, Nm=N): return integrate(eps_integrand(g, Nm)).real

if __name__ == '__main__':
    import sys, time, numpy as np
    sys.path.insert(0, '.')
    print("C.15 Arb callbacks self-test (NUMERICAL comparison against float pipeline; prec =", ctx.prec, "bits)")
    print("q_T =", QT.str(20), " (c15_entries_ref: 1.27619777770742)")
    ref = [[0.207950766678, -0.381518175204, 0.220969487769, -0.0975943904542],
           [-0.381518175204, 0.892990625409, -0.14157757611, 0.0410416374567],
           [0.220969487769, -0.14157757611, 1.20151232126, 0.0260117645367],
           [-0.0975943904542, 0.0410416374567, 0.0260117645367, 1.39953509124]]
    t0 = time.time(); mx = 0
    for (j, k) in [(0, 0), (0, 1), (1, 2), (2, 3), (3, 3)]:
        v = A_entry(j, k); d = abs(float(v.mid()) - ref[j][k]); mx = max(mx, d)
        print(f"A[{j},{k}] = {v.str(18)}   ref {ref[j][k]}   |diff| = {d:.2e}")
    print("max |A_arb - A_ref(mpmath 12 digits)| =", f"{mx:.2e}", f"({time.time()-t0:.1f}s)")
    from c16_schur_general import band_matrix, mabs, weighted_tail
    Lf = float(L.mid()); P = [(2, float(C2.mid()))]
    Mf = band_matrix(8, Lf, 60.0, 'even', lambda t: mabs(t, 60.0, P))
    for (j, k) in [(0, 0), (0, 5), (7, 7)]:
        v = M_entry(j, k); print(f"M[{j},{k}] = {v.str(18)}   float {Mf[j,k]:.12f}   |diff| = {abs(float(v.mid())-Mf[j,k]):.2e}")
    t0 = time.time()
    for name, g, gf in [("m_abs", m_abs, lambda t: mabs(t, 60.0, P)), ("m_neg", m_neg, lambda t: 0*t + float(OMT.mid()) - __import__('c16_schur_general').Om(t))]:
        v = eps_entry(g); f = weighted_tail(N, Lf, 60.0, 'even', gf)
        print(f"eps^{name}_64 (closed-form K_W, Arb) = {v.str(15)}   float mode-sum+remainder {f:.6e}   rel diff {abs(float(v.mid())-f)/f:.1e}")
    print(f"({time.time()-t0:.1f}s)")
    print("sharp |a_T| callback with analytic=True returns:", sharp_abs_integrand(0, 0)(acb(1), True))
