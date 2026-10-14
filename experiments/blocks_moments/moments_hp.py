# High-precision A_k = (1/2pi) int_R P(l) D(l)^k dl, P=pi l/sinh(pi l), D=Re psi(1/2+il/2)-psi(1/2)
# float mpmath quadrature (NOT interval-certified); precision checked by two working precisions.
from mpmath import mp, mpf, pi, sinh, digamma, re, quad, log, zeta, polylog, pslq, nstr, inf, mpc
def A(k, dps):
    mp.dps = dps
    P = lambda l: pi*l/sinh(pi*l) if l != 0 else mpf(1)
    D = lambda l: re(digamma(mpf(1)/2 + mpc(0,l)/2)) - digamma(mpf(1)/2)
    return quad(lambda l: P(l)*D(l)**k, [0,1,2,4,8,16,32,64,inf])/pi
res = {}
for k in range(0,6):
    a = A(k, 75); b = A(k, 90); mp.dps = 75
    res[k] = b
    print(f"A_{k} = {nstr(b,60)}   |diff 75 vs 90 dps| = {nstr(abs(a-b),3)}")
mp.dps = 70
L = log(2)
print("A_3 - (log2/2 - 3/16) =", nstr(res[3] - (L/2 - mpf(3)/16), 5))
Li4 = polylog(4, mpf(1)/2)
names = ["1","log2","log2^2","pi^2","zeta3","log2^3","pi^2log2","pi^4","zeta3log2","pi^2log2^2","log2^4","Li4(1/2)"]
basis = [mpf(1),L,L**2,pi**2,zeta(3),L**3,pi**2*L,pi**4,zeta(3)*L,pi**2*L**2,L**4,Li4]
for k in (3,4,5):
    mp.dps = 60
    for nb in (5, 8, 12):
        r = pslq([res[k]] + basis[:nb], maxcoeff=10**6, maxsteps=10**6)
        print(f"PSLQ A_{k} basis[:{nb}]:", r and list(zip(["A"]+names[:nb], r)))
