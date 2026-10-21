"""Grok S70 review: independent check of Codex constants (Oct 20/21). NUMERICAL (mpmath), not interval."""
import mpmath as mp
mp.mp.dps=40
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
x0=xi(mp.mpf(1)/2); x2=mp.diff(xi,mp.mpf(1)/2,2)
r=x2/x0; print("xi''(1/2)/xi(1/2) =",mp.nstr(r,25),"(Codex 0.0462099862308379415778676...)")
d1=mp.diff(lambda s: mp.log(xi(s)),1); print("xi'(1)/xi(1) =",mp.nstr(d1,25)," closed 1+g/2-log2-log(pi)/2 =",mp.nstr(1+mp.euler/2-mp.log(2)-mp.log(mp.pi)/2,25))
Cxi=8*(r/2-d1); print("C_xi =",mp.nstr(Cxi,25),"(Codex 0.0000742731943834957969885...)")
# zero-sum cross-check with first 300 zeros (partial)
s=mp.mpf(0); s2=mp.mpf(0)
for k in range(1,301):
    g=mp.im(mp.zetazero(k)); s+=2/(g**2*(g**2+mp.mpf(1)/4)); s2+=2/g**2
print("300-zero partial 2sum 1/(g^2(g^2+1/4)) =",mp.nstr(s,15)," 2sum 1/g^2 =",mp.nstr(s2,15))
kap=2+mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi); print("kappa =",mp.nstr(kap,20))
# C_star = -2gamma -2U -2 + H_inf, U=sum_p log p/(p(p-1)), H_inf=2 sum_p log p p^{-3/2}/(1-p^{-1/2})
from sympy import primerange
P=10**7; U=mp.mpf(0); H=mp.mpf(0)
for p in primerange(2,P):
    lp=mp.log(p); U+=lp/(p*(p-1)); H+=2*lp*mp.mpf(p)**(-1.5)/(1-mp.mpf(p)**(-0.5))
# tails via PNT density: sum_{p>P} log p f(p) ~ int_P^inf f(x)dx
Ut=1/mp.mpf(P); Ht=2*2*mp.mpf(P)**(-0.5)
print("U(1e7)=",mp.nstr(U,15),"+tail~",mp.nstr(Ut,3)," H_inf(1e7)=",mp.nstr(H,15),"+tail~",mp.nstr(Ht,3))
U2=U+Ut; H2=H+Ht; print("C_star ~",mp.nstr(-2*mp.euler-2*U2-2+H2,12),"(Codex partial-1e6 0.2845084973; tail-corrected est. above)")
print("known: -gamma-U = E =",mp.nstr(-mp.euler-U2,12),"(lit. -1.3325822757)")
