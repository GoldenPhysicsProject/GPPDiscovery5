# Checks: (1) P = Kubo-Mori kernel at beta=2pi; (2) archimedean Weil term as Levy-Khintchine
# with jump density nu(u)=2e^{-u/2}/(1-e^{-2u}) = e^{u/2}/sinh(u); (3) explicit formula with
# arch+prime as ONE positive jump measure, vs zeros. Float (mpmath 30 dps), not certified.
from mpmath import mp,mpf,pi,sinh,exp,cos,psi,quad,inf,zetazero,log,euler,mangoldt,sqrt
mp.dps=30
# (1)
for lam in [mpf('0.7'),mpf('2.3')]:
    kap=(2*pi*lam/2)/sinh(2*pi*lam/2); print('KM kernel - P', kap-pi*lam/sinh(pi*lam))
# (2) Re psi(1/4+ir/2)-psi(1/4) = int_0^inf (1-cos ru) nu(u) du
nu=lambda u: 2*exp(-u/2)/(1-exp(-2*u))
for r in [mpf(1),mpf(5),mpf(14)]:
    lhs=psi(0,mpf(1)/4+1j*r/2).real-psi(0,mpf(1)/4)
    rhs=quad(lambda u:(1-cos(r*u))*nu(u),[0,1,5,20,inf],maxdegree=10)
    print('Levy r=',r,lhs-rhs)
# (3) g(u)=exp(-u^2/(2s^2)), H(r)=s sqrt(2pi) exp(-s^2 r^2/2)
s=mpf('0.5'); g=lambda u: exp(-u**2/(2*s*s)); H=lambda r: s*sqrt(2*pi)*exp(-s*s*r*r/2)
Z=2*sum(H(zetazero(k).imag) for k in range(1,40))
pole=quad(lambda u:g(u)*(exp(u/2)+exp(-u/2)),[-inf,0,inf])
jump=quad(lambda u:(exp(-2*u)*g(0)*2/ (2*exp(-u/2)) *nu(u)*0+ (exp(-2*u)*g(0)-exp(-u/2)*g(u))*2/(1-exp(-2*u))),[0,1,inf])
prime=2*sum(mangoldt(n)*n**mpf(-0.5)*g(log(n)) for n in range(2,200))
geo=pole-(euler+log(pi))*g(0)+jump-prime
print('zeros',Z,'geometric',geo,'diff',Z-geo)
print('arch jump mass total vs prime mass in window: arch part',jump,'prime',prime)
