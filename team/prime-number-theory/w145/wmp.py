from mpmath import mp, mpf, exp, log, euler, zeta, sqrt, diff, quad
import math
mp.dps=40
X=10**5
lam={}
for n in range(2,X+1):
    pass
# prime powers
sv=[True]*(X+1); sv[0]=sv[1]=False
for p in range(2,int(X**.5)+1):
    if sv[p]:
        for k in range(p*p,X+1,p): sv[k]=False
pp=[]
for p in range(2,X+1):
    if sv[p]:
        q=p
        while q<=X: pp.append((q,log(mpf(p)))); q*=p
pp.sort()
g=+euler
def F(t,A,B):  # antiderivative of r^2, r=2e^{t/2}(A-t)-B
    u=A-t
    return 4*exp(t)*(u*u+2*u+2) -4*B*(2*exp(t/2)*(u+2)) + B*B*t
tot=mpf(0); S1=mpf(0); S2=mpf(0); a=mpf(0)
for (n,L) in pp+[(X,None)]:
    b=log(mpf(n))
    A=S2+2+g; B=2*S1+2
    tot+=F(b,A,B)-F(a,A,B)
    if L is not None:
        S1+=L/sqrt(n); S2+=L/n
    a=b
print("W(log 1e5) mpmath dps40 closed-form:",mp.nstr(tot,15))
Jinf=2-2*diff(zeta,mpf(1)/2)/zeta(mpf(1)/2)
print("J_inf=2-2 zeta'/zeta(1/2)=",mp.nstr(Jinf,20)," J_inf^2=",mp.nstr(Jinf**2,15))
print("-2 zeta'(1/2)=",mp.nstr(-2*diff(zeta,mpf(1)/2),30))
