"""C.16/C.17 Arb callbacks, general L, prime-power set S and parity (odd p=2 certificate of Thm C.16.A; reusable by Prime).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN.
Same conventions and domain/analyticity choices (1)-(6) as c15_arb_callbacks.py; differences:
  * odd modes: e_k = L^{-1/2} sin(w_k x), w_k = k pi/L (k >= 1; index k=0..N-1 below means w = (k+1) pi/L),
    real transform ehat_k = sqrt(L)[sinc((w-t)L) - sinc((w+t)L)] (entire);
  * odd full kernel (Lemma C.16.K): K_o(t) = [L - L sinc(2tL) - 2 sin(tL)(coth(L/2) sin(tL)/2 - t cos(tL))/(t^2+1/4)]/(1/4+t^2),
    meromorphic (poles at +-i/2 only); cancels near t=0 (K_o(0)=0) and against the N-mode sum: use prec >= 128;
  * q(t) = Omega(t) - sum_{n in S} 2 c_n cos(t log n), c_n = Lambda(n)/sqrt(n); q_T = Omega(T) - 2 sum c_n;
    m_abs = Omega(T) - Omega(t) + sum 2 c_n (1 - cos(t log n)), m_neg = Omega(T) - Omega(t)  (all meromorphic).
Usage: cb = Callbacks(twoL=arb(3).log(), S=[(2, arb(2).log()/arb(2).sqrt())], parity='odd', T=60, N=128); cb.A_entry(j,k) ...
"""
from flint import acb, arb, ctx
ctx.prec = 128
PI = arb.pi(); LOGPI = PI.log(); Q4 = arb(1)/4
def Omega(t):
    h = acb(0, 1)*t/2
    return ((Q4 + h).digamma() + (Q4 - h).digamma())/2 - LOGPI
class Callbacks:
    def __init__(s, twoL, S, parity, T, N, panel=5):
        s.L = arb(twoL)/2; s.S = [(arb(n).log(), arb(c)) for n, c in S]; s.par = parity; s.T = arb(T); s.N = N; s.panel = panel
        s.OMT = Omega(acb(s.T)).real; s.QT = s.OMT - 2*sum(c for _, c in s.S)
        s.sg = 1 if parity == 'even' else -1
    def w(s, k): return ((arb(k) + arb(1)/2) if s.par == 'even' else arb(k + 1))*PI/s.L
    def nk(s, k): return (s.w(k)**2 + Q4).sqrt()
    def ehat(s, k, t):
        wk = s.w(k); return s.L.sqrt()*(((wk - t)*s.L).sinc() + s.sg*((wk + t)*s.L).sinc())
    def q(s, t): return Omega(t) - sum(2*c*(lg*t).cos() for lg, c in s.S)
    def m_abs(s, t): return s.OMT - Omega(t) + sum(2*c*(1 - (lg*t).cos()) for lg, c in s.S)
    def m_neg(s, t): return s.OMT - Omega(t)
    def Kfull(s, t):
        L = s.L; c, sn = (t*L).cos(), (t*L).sin()
        if s.par == 'even': num = L + L*(2*t*L).sinc() - 2*c*((L/2).tanh()/2*c + t*sn)/(t*t + Q4)
        else: num = L - L*(2*t*L).sinc() - 2*sn*(1/(L/2).tanh()/2*sn - t*c)/(t*t + Q4)
        return num/(Q4 + t*t)
    def KW(s, t): return s.Kfull(t) - sum(s.ehat(k, t)**2/(s.w(k)**2 + Q4) for k in range(s.N))
    def A_integrand(s, j, k):
        f = 1/(PI*s.nk(j)*s.nk(k)); return lambda t, an: f*(Q4 + t*t)*(s.q(t) - s.QT)*s.ehat(j, t)*s.ehat(k, t)
    def M_integrand(s, j, k):
        f = 1/(PI*s.nk(j)*s.nk(k)); return lambda t, an: f*(Q4 + t*t)*s.m_abs(t)*s.ehat(j, t)*s.ehat(k, t)
    def eps_integrand(s, g): return lambda t, an: (Q4 + t*t)*g(t)*s.KW(t)/PI
    def integrate(s, f):
        tot = acb(0); x = arb(0)
        while x < s.T:
            y = x + s.panel
            if y > s.T: y = s.T
            tot += acb.integral(f, x, y); x = y
        return tot.real
    def A_entry(s, j, k): return s.integrate(s.A_integrand(j, k)) + (s.QT if j == k else 0)
    def M_entry(s, j, k): return s.integrate(s.M_integrand(j, k))
    def eps(s, which): return s.integrate(s.eps_integrand(s.m_abs if which == 'abs' else s.m_neg))
if __name__ == '__main__':
    import sys, numpy as np; sys.path.insert(0, '.')
    from c16_schur_general import band_matrix, qsym, qT_of, mabs, mneg, weighted_tail, prime_set
    c2 = arb(2).log()/arb(2).sqrt()
    cb = Callbacks(arb(3).log(), [(2, c2)], 'odd', 60, 128)
    Lf = float(cb.L.mid()); P = prime_set(Lf); qT = qT_of(60., P)
    Af = band_matrix(8, Lf, 60., 'odd', lambda t: qsym(t, P) - qT) + qT*np.eye(8)
    Mf = band_matrix(8, Lf, 60., 'odd', lambda t: mabs(t, 60., P))
    print("odd, 2L = log 3, T = 60, N = 128 (self-test, NUMERICAL):  q_T =", cb.QT.str(18))
    for j, k in [(0, 0), (0, 1), (2, 5), (7, 7)]:
        a = cb.A_entry(j, k); m = cb.M_entry(j, k)
        print(f"A[{j},{k}] = {a.str(16)} float {Af[j,k]:.13f} diff {abs(float(a.mid())-Af[j,k]):.1e} | M[{j},{k}] = {m.str(16)} float {Mf[j,k]:.13f} diff {abs(float(m.mid())-Mf[j,k]):.1e}", flush=True)
    for wh, g in (('abs', lambda t: mabs(t, 60., P)), ('neg', lambda t: mneg(t, 60., P))):
        e = cb.eps(wh); f = weighted_tail(128, Lf, 60., 'odd', g)
        print(f"eps^{wh}_128 Arb (closed-form K_o) = {e.str(12)}  float mode-sum+remainder (upper bound) {f:.6e}", flush=True)
