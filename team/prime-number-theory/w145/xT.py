import numpy as np, math
from mpmath import euler
g=float(euler); N=2*10**6
lam=np.zeros(N+1); sv=np.ones(N+1,bool); sv[:2]=False
for p in range(2,int(N**.5)+1):
    if sv[p]: sv[p*p::p]=False
for p in np.nonzero(sv)[0]:
    q=p
    while q<=N: lam[q]=math.log(p); q*=p
pp=np.nonzero(lam)[0]; S1=np.cumsum(lam[pp]/np.sqrt(pp)); S2=np.cumsum(lam[pp]/pp)
xs,ws=np.polynomial.legendre.leggauss(8)
def stats(T):
    X=math.exp(T); m=pp<=X
    e=np.concatenate(([1.0],pp[m].astype(float),[X])); s1=np.concatenate(([0.],S1[m])); s2=np.concatenate(([0.],S2[m]))
    a=np.log(e[:-1]); b=np.log(e[1:]); W=0; x=0
    for xi,wi in zip(xs,ws):
        t=(a+b)/2+(b-a)/2*xi; r=2*np.exp(t/2)*(s2-t+2+g)-2*s1-2; w=wi*(b-a)/2
        W+=np.sum(w*r*r); x+=np.sum(w*np.exp(-(T-t)/2)*r)
    return W,x
for T in [5,7.27,8,11,14.5,math.log(2e6)]:
    W,x=stats(T); print(f"T={T:.4f} W={W:.4f} x(T)={x:.5f} x^2={x*x:.3f} E_c=W-x^2={W-x*x:.3f}")
