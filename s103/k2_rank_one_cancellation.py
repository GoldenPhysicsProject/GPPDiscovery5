"""S103 (Grok): test of ARP G.3d. The k=2 rank-one mean (+1/2 * int h) is cancelled
exactly in the mean by the -sqrt(x) secondary term of theta(x) inside the k=1 layer.
Numerics only test this one derivation step (Daniel, Oct 22 pivot). RH OPEN."""
import numpy as np
from mpmath import mp, quad, exp
U=16.0                      # h supported on [0,U]; n <= e^U ~ 8.9e6
N=int(np.exp(U))+1
s=np.ones(N+1,bool); s[:2]=False
for i in range(2,int(N**.5)+1):
    if s[i]: s[i*i::i]=False
p=np.nonzero(s)[0].astype(float); lp=np.log(p)
def h(u):   # smooth bump on (2,U)
    u=np.asarray(u,float); a,b=2.0,U; out=np.zeros_like(u)
    m=(u>a)&(u<b); x=(u[m]-a)/(b-a); out[m]=np.exp(-1/(x*(1-x)))*np.exp(4); return out
def hm(u):
    a,b=2.0,U
    if not (a<u<b): return 0
    x=(u-a)/(b-a); return exp(-1/(x*(1-x)))*exp(4)
mp.dps=30
I0=float(quad(hm,[2,9,U]))                         # int h
Ih=float(quad(lambda u: hm(u)*exp(u/2),[2,9,U]))   # PNT mean of k=1 (pole partner)
S1=np.sum(lp*p**-0.5*h(lp))      # k=1 layer
S2=np.sum(lp/p*h(2*lp))          # k=2 layer
pk=[np.sum(lp*p**(-k/2)*h(k*lp)) for k in range(3,25)]
print("int h =",I0)
print("k=1 fluctuation  S1 - int h e^{u/2} =",S1-Ih, " (prediction: -1/2 int h =",-I0/2,"+ zero oscillation)")
print("k=2 fluctuation  S2 - (1/2) int h  =",S2-I0/2)
print("k=2 layer S2 =",S2," vs 1/2 int h =",I0/2)
print("Lambda total minus PNT mean (all k):",S1+S2+sum(pk)-Ih)
print("k>=3 total:",sum(pk))
