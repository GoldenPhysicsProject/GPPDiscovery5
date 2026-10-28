# N.13 single-step constant evaluation (NUMERICAL check of closed forms; not a proof step).
from mpmath import mp, euler, pi, log, exp, sqrt, digamma, mpf, zeta
mp.dps = 20
c2 = log(2)/sqrt(2)
Om0 = -euler - pi/2 - 3*log(2) - log(pi)
print("c2", c2, "Omega(0) closed form", Om0, "check", digamma(mpf(1)/4) - log(pi))
t = mpf('12.8')
lb = Om0 + 4*t**2/(mpf(1)/4 + t**2) + mpf('0.8')*t**2/(mpf('6.25') + t**2) + log(1 + t**2/mpf('20.25'))/2
print("lower bound Omega(12.8)", lb, ">= c2:", lb >= c2)
print("zeta(3,1/4) <= 64.84:", zeta(3, mpf(1)/4))
y = mpf('0.186'); print("h* crossing", 1 + mpf('0.6495')/y, 2*mpf('64.84')*y**2)
X0 = abs(Om0) + 2*c2
print("X(m0) = %s + %s/m0" % (X0, 4*c2*abs(Om0)))
print("ell coefficient pi/(16 c2 e^X0) =", pi/(16*c2*exp(X0)))
print("dilation slope (2 + 4.5 + 2|Om0|)/log2 =", (2 + mpf('4.5') + 2*abs(Om0))/log(2))
