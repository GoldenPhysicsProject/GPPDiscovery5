# NUMERICAL (tiny step check for §N(c)): Omega(0) = Re psi(1/4) - log pi and m_X(0) = Omega(0) - 2*sum_{p<=X} log p/(sqrt p - 1).
# Also Astra's normalisation constants int b = 1/2, int G0 = (int b)^2 = 1/4 for b = sin^2(pi u) on [0,1].
from mpmath import mp, digamma, log, pi, sqrt, quad, sin, euler
def primerange(a,b):
    sv=bytearray([1])*b; sv[0:2]=b'\x00\x00'
    for i in range(2,int(b**0.5)+1):
        if sv[i]: sv[i*i::i]=bytearray(len(sv[i*i::i]))
    return [i for i in range(a,b) if sv[i]]
mp.dps = 30
Om0 = digamma(mp.mpf(1)/4).real - log(pi)
print("psi(1/4) =", digamma(mp.mpf(1)/4), " closed form -gamma-pi/2-3log2 =", -euler - pi/2 - 3*log(2))
print("Omega(0) =", Om0)
for X in (2, 23, 10**3, 10**5):
    A = 2*sum(log(p)/(sqrt(p)-1) for p in primerange(2, X+1))
    print("X =", X, " 2*D_X(0) =", A, " m_X(0) =", Om0 - A)
ib = quad(lambda u: sin(pi*u)**2, [0, 1])
print("int b =", ib, " (int b)^2 =", ib**2)
