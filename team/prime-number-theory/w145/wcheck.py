import numpy as np, math
from mpmath import mp, mpf, euler
N=2*10**6
# von Mangoldt
lam=np.zeros(N+1)
sieve=np.ones(N+1,bool); sieve[:2]=False
for p in range(2,int(N**.5)+1):
    if sieve[p]: sieve[p*p::p]=False
primes=np.nonzero(sieve)[0]
for p in primes:
    lp=math.log(p); q=p
    while q<=N: lam[q]=lp; q*=p
pp=np.nonzero(lam)[0]            # prime powers (cell boundaries)
g=float(euler)
S1=np.cumsum(lam[pp]/np.sqrt(pp)); S2=np.cumsum(lam[pp]/pp); PSI=np.cumsum(lam[pp])
# cells: [1,pp0) with S=0, then [pp_i, pp_{i+1})
def W_of(X, tau=0.0, npts=6):
    xs,ws=np.polynomial.legendre.leggauss(npts)
    edges=np.concatenate(([1.0],pp[pp<=X].astype(float),[X]))
    s1=np.concatenate(([0.0],S1[pp<=X])); s2=np.concatenate(([0.0],S2[pp<=X]))
    a=np.log(edges[:-1]); b=np.log(edges[1:])
    tot=0.0
    for xi,wi in zip(xs,ws):
        t=(a+b)/2+(b-a)/2*xi
        r=2*np.exp(t/2)*(s2-t+2+g+tau)-2*s1-2
        tot+=np.sum(wi*(b-a)/2*r*r)
    return tot
X=2e6; T=math.log(X)
print("T=log 2e6",T)
print("W exact-tail, X=2e6:",W_of(X))
print("W exact-tail, T=14.5:",W_of(math.exp(14.5)))
# truncated tail variant: future integral of q cut at Tmax=log N
s2N=S2[-1]; psiN=PSI[-1]
tau=(-1-g)-(s2N-psiN/N-math.log(N))   # = int_{Tmax}^inf q
print("tau=int_Tmax^inf q =",tau)
# dropping it means replacing int_t^inf q by int_t^inf q - tau -> r changes by +2 sqrt(x) tau
print("W tail-dropped, X=2e6:",W_of(X,tau=+tau))
print("W tail-dropped, T=14.5:",W_of(math.exp(14.5),tau=+tau))
print("W tail-dropped, T=7.27:",W_of(math.exp(7.27),tau=+tau), " exact:",W_of(math.exp(7.27)))
# boundary term x(T): x' + x/2 = r. E_c = W - x^2
