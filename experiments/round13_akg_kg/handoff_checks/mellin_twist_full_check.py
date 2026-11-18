"""Heart arm: independent check of the FULL balanced Mellin twist (D2 Result E) on non-constant tests.
W_theta(f) == W_zeta(cosh(theta x) f) - W_zeta(sinh(theta x) f), zero-free source forms only.
Normalisation (real f supported in [-L/2,L/2], A_f(t)=int f(x)f(x+t)dx, w(t)=e^{-t/2}/(1-e^{-2t})):
 W_theta(f) = C_{inf,theta} A(0) + 2 int_0^inf w(t)cosh(theta t)[A(0)-A(t)]dt
            + 2 int A(t) cosh(t/2) cosh(theta t) dt - 2 sum_n Lambda(n) n^{-1/2} cosh(theta log n) A(log n),
 C_{inf,theta} = (psi(1/4+theta/2)+psi(1/4-theta/2))/2 - log pi;  W_zeta = W_0.
No zeta zeros are used anywhere.
"""
import mpmath as mp, json
mp.mp.dps = 20
def vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return mp.log(p) if m == 1 else mp.mpf(0)
    return mp.mpf(0)
def W(fun, L, th):
    th = mp.mpf(th)
    A = lambda t: mp.quad(lambda x: fun(x)*fun(x+t), [-L/2, L/2-t]) if t < L else mp.mpf(0)
    A0 = A(0)
    C = (mp.digamma(mp.mpf(1)/4+th/2)+mp.digamma(mp.mpf(1)/4-th/2))/2 - mp.log(mp.pi)
    w = lambda t: mp.exp(-t/2)/(-mp.expm1(-2*t))
    arch = C*A0 + 2*mp.quad(lambda t: w(t)*mp.cosh(th*t)*(A0-A(t)), [0, L/4, L/2, L]) \
           + 2*A0*mp.quad(lambda t: w(t)*mp.cosh(th*t), [L, mp.inf])
    pole = 2*2*mp.quad(lambda t: A(t)*mp.cosh(t/2)*mp.cosh(th*t), [0, L/2, L])
    prime = 2*mp.fsum(vm(n)/mp.sqrt(n)*mp.cosh(th*mp.log(n))*A(mp.log(n)) for n in range(2, int(mp.exp(L))+1))
    return arch + pole - prime
out = []
for L, th in [(1.5, 0.2), (2.0, 0.35), (2.5, 0.45)]:
    L = mp.mpf(L)
    f = lambda x, L=L: (1-(2*x/L)**2)**2*(1+mp.mpf('0.3')*x+mp.mpf('0.2')*x**2) if abs(x) <= L/2 else mp.mpf(0)
    lhs = W(f, L, th)
    rhs = W(lambda x: mp.cosh(th*x)*f(x), L, 0) - W(lambda x: mp.sinh(th*x)*f(x), L, 0)
    out.append(dict(L=float(L), theta=th, W_theta=mp.nstr(lhs, 20), rhs=mp.nstr(rhs, 20), absdiff=mp.nstr(abs(lhs-rhs), 3)))
    print(out[-1])
json.dump(out, open('mellin_twist_full_check.json', 'w'), indent=1)
