"""W(T) = int_0^T e^t (j(t) + 2(1+gamma))^2 dt, exact cell-by-cell (mpmath).
Convention (Astra / Grok S77): T = log X; psi = full von Mangoldt sum (prime powers included);
j = J_{1/2} with the exact future tail, so that the integrand at time t is
    r(t) = e^{t/2}(j(t) + 2(1+gamma)) = 2 e^{t/2} (S2(t) - t + 2 + gamma) - 2 S1(t) - 2,
    S1(t) = sum_{n<=e^t} Lambda(n) n^{-1/2},  S2(t) = sum_{n<=e^t} Lambda(n)/n.
r is piecewise of the form 2e^{t/2}(A - t) - B between consecutive prime powers, so W, x(T) and E_c
have closed-form antiderivatives per cell: no quadrature, no truncated tail.
Also x(T) = e^{-T/2} int_0^T e^{s/2} r(s) ds (x' + x/2 = r, x(0)=0) and E_c(T) = W(T) - x(T)^2.
Usage: python w_exact.py T [dps]     e.g.  python w_exact.py 14.5 ;  python w_exact.py 'log(2e6)'
NUMERICAL (mpmath, working precision dps; the Arb-CERTIFIED values are Grok's s77astra/zero_damping_check.py).
"""
import sys, math
from mpmath import mp, mpf, exp, log, sqrt, euler
mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 30
Tin = sys.argv[1]
T = log(mpf(2)*10**6) if Tin == 'log(2e6)' else (log(mpf(10)**5) if Tin == 'log(1e5)' else mpf(Tin))
X = int(math.floor(float(exp(T)))) + 1
sv = bytearray([1]) * (X + 1); sv[0] = sv[1] = 0
for p in range(2, int(X ** .5) + 1):
    if sv[p]: sv[p*p::p] = bytearray(len(range(p*p, X + 1, p)))
pp = []
for p in range(2, X + 1):
    if sv[p]:
        q = p
        while q <= X: pp.append((q, p)); q *= p
pp.sort()
g = +euler
def FW(t, A, B):          # antiderivative of r^2
    u = A - t
    return 4*exp(t)*(u*u + 2*u + 2) - 8*B*exp(t/2)*(u + 2) + B*B*t
def FX(t, A, B):          # antiderivative of e^{t/2} r
    return 2*exp(t)*(A - t + 1) - 2*B*exp(t/2)
W = mpf(0); I = mpf(0); S1 = mpf(0); S2 = mpf(0); a = mpf(0)
for (n, p) in pp:
    b = log(mpf(n))
    if b > T: break
    A = S2 + 2 + g; B = 2*S1 + 2
    W += FW(b, A, B) - FW(a, A, B); I += FX(b, A, B) - FX(a, A, B)
    lp = log(mpf(p)); S1 += lp/sqrt(n); S2 += lp/n; a = b
A = S2 + 2 + g; B = 2*S1 + 2
W += FW(T, A, B) - FW(a, A, B); I += FX(T, A, B) - FX(a, A, B)
x = exp(-T/2)*I
print("T=%s  W=%s  x(T)=%s  E_c=%s" % (mp.nstr(T, 12), mp.nstr(W, 15), mp.nstr(x, 12), mp.nstr(W - x*x, 12)))
