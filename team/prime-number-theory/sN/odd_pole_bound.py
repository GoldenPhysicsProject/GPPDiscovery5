# NUMERICAL evaluation of the constants in Prop. N.3 (the inequality is PROVED; this only evaluates it).
# f_c(u) = b(u-c) - b(u+c), b(u) = cos^2(pi u/(2w)) on [-w,w]. Q^{P_X}(f_c) < 0 whenever
#   8 sinh^2(c/2) B^2 > 4 (A_X ||b||^2 + omega_b),  B = int b e^{u/2},
#   omega_b = (1/2pi) int |bhat|^2 Omega^+ dt,  A_X = 2 sum_{p<=X} log p/(sqrt p - 1),
#   Omega(t) = Re digamma(1/4+it/2) - log pi.
from mpmath import mp, mpf, quad, sin, pi, exp, log, sqrt, digamma, asinh, findroot, sinh
mp.dps = 15
def primes(n):
    s = bytearray([1])*(n+1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
def A(X): return 2*sum(log(p)/(sqrt(p)-1) for p in primes(X))
def Om(t): return digamma(mpf(1)/4 + 0.5j*t).real - log(pi)
t0 = findroot(Om, 8)
def S(a, w): return 2*w if a == 0 else 2*sin(a*w)/a        # int_{-w}^{w} cos(a u) du
def bhat(t, w):  # b = (1 + cos(pi u/w))/2 on [-w,w]
    k = pi/w
    return S(t, w)/2 + (S(t-k, w) + S(t+k, w))/4
def consts(w):
    B = quad(lambda u: (1+mp.cos(pi*u/w))/2*exp(u/2), [-w, w])
    nb = quad(lambda u: ((1+mp.cos(pi*u/w))/2)**2, [-w, w])
    om = 2*quad(lambda t: bhat(t, w)**2*Om(t), [t0, 50, 200, 2000])/(2*pi)
    return B, nb, om
print("Omega zero t0 =", t0)
cs = {w: consts(mpf(w)) for w in (0.75, 1.0, 1.25, 1.5, 2.0)}
for X in (23, 10**3, 10**5, 10**7):
    AX = A(X); best = None
    for w, (B, nb, om) in cs.items():
        c = 2*asinh(sqrt((AX*nb + om)/(2*B**2)))
        if best is None or c + w < best[0]: best = (c + w, w, c, B, nb, om)
    L, w, c, B, nb, om = best
    print(f"X={X} A_X={float(AX):.4f} w={w} c*={float(c):.4f} L1_bound={float(L):.4f} (log X)/2={float(log(X)/2):.4f} "
          f"excess={float(L - log(X)/2):.4f} B={float(B):.4f} |b|^2={float(nb):.4f} omega_b={float(om):.3e}")
