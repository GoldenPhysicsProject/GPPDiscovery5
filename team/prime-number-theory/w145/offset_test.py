"""Diagnosis by elimination: does a constant offset r -> r - delta (or a scale r -> s*r) reproduce ALL of Turbo's
reported numbers (W(14.5)=158.1, W-x^2=114.20, x(T)^2 = 36.7, 42.1, 43.4, 43.9 at T=5,8,11,14.5)? NUMERICAL (float64, fine grid)."""
import numpy as np, math
g=0.5772156649015329; N=2*10**6
lam=np.zeros(N+1); sv=np.ones(N+1,bool); sv[:2]=False
for p in range(2,int(N**.5)+1):
    if sv[p]: sv[p*p::p]=False
for p in np.nonzero(sv)[0]:
    q=p
    while q<=N: lam[q]=math.log(p); q*=p
n=np.arange(N+1); S1c=np.cumsum(lam/np.sqrt(np.maximum(n,1))); S2c=np.cumsum(lam/np.maximum(n,1))
h=2e-5; t=np.arange(0,14.5+h/2,h); idx=np.floor(np.exp(t)+1e-12).astype(int)
r=2*np.exp(t/2)*(S2c[idx]-t+2+g)-2*S1c[idx]-2
def stats(rr):
    W=np.cumsum(rr**2)*h; I=np.cumsum(np.exp(t/2)*rr)*h; x=np.exp(-t/2)*I
    out={}
    for T in (5,8,11,14.5):
        k=int(round(T/h)); out[T]=(W[k],x[k]**2)
    return out
turbo={'W':158.1,'bulk':114.20,'x2':[36.7,42.1,43.4,43.9]}
base=stats(r)
print("exact  : W=%.3f bulk=%.3f x2=%s"%(base[14.5][0],base[14.5][0]-base[14.5][1],[float(round(base[T][1],2)) for T in (5,8,11,14.5)]))
print("ratios Turbo/exact: W %.4f bulk %.4f x2 %s"%(158.1/base[14.5][0],114.2/(base[14.5][0]-base[14.5][1]),[float(round(a/base[T][1],4)) for a,T in zip(turbo['x2'],(5,8,11,14.5))]))
for d in (0.04,0.05,0.06,0.07):
    s=stats(r-d); print("offset %.2f: W=%.2f bulk=%.2f x2=%s"%(d,s[14.5][0],s[14.5][0]-s[14.5][1],[float(round(s[T][1],2)) for T in (5,8,11,14.5)]))
for sc in (0.980,0.981,0.982):
    s=stats(sc*r); print("scale %.3f: W=%.2f bulk=%.2f x2=%s"%(sc,s[14.5][0],s[14.5][0]-s[14.5][1],[float(round(s[T][1],2)) for T in (5,8,11,14.5)]))
