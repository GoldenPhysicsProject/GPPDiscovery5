"""Verify Celestial's zero of Z_Q(s)=sum' Q(m,n)^{-s}, Q=2x^2+2xy+3y^2 (disc -20), via
Z_Q = zeta*L(chi_-20) - L(chi_-4)*L(chi_5)   (r_1+r_2 = 2 sum chi_-20(d), r_1-r_2 = 2 sum chi_-4(d)chi_5(n/d); w=2).
Standard normalization: converges Re s>1, Lambda(s)=(sqrt20/2pi)^s Gamma(s) Z(s) = Lambda(1-s), critical line 1/2."""
import mpmath as mp
mp.mp.dps=30
def kron(D,n):
    return int(mp.mpf(0)) if False else None
chi20=[0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0]  # placeholder, computed below
def kr(D,n):
    # Kronecker symbol (D/n) for n>0 via sympy-free implementation
    from math import gcd
    if gcd(D,n)!=1: return 0
    res=1; a=D; m=n
    # factor out 2 from m
    while m%2==0:
        m//=2
        if a%8 in (3,5): res=-res
    # Jacobi (a/m) with a possibly negative
    if a<0:
        a=-a
        if m%4==3: res=-res
    a%=m
    while a:
        while a%2==0:
            a//=2
            if m%8 in (3,5): res=-res
        a,m=m,a
        if a%4==3 and m%4==3: res=-res
        a%=m
    return res if m==1 else 0
c20=[kr(-20,n) for n in range(20)]; c4=[kr(-4,n) for n in range(4)]; c5=[kr(5,n) for n in range(5)]
Zq=lambda s: mp.zeta(s)*mp.dirichlet(s,c20)-mp.dirichlet(s,c4)*mp.dirichlet(s,c5)
# sanity: direct sum at s=3
S=mp.mpf(0)
for x in range(-60,61):
    for y in range(-60,61):
        if x or y: S+=mp.mpf(2*x*x+2*x*y+3*y*y)**-3
print("direct s=3",S,"formula",Zq(3))
for s0 in [mp.mpc(2.079,7.629), mp.mpc(2.079/2,7.629/2), mp.mpc(2.079/2+0.5,7.629/2)]:
    print(s0, abs(Zq(s0)))
r=mp.findroot(Zq, mp.mpc(2.079,7.629)); print("root",r, abs(Zq(r)))
# also check is it a zero of the x^2+5y^2 one or zeta_K
