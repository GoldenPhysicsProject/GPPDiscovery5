# NUMERICAL: semi-local truncated Weil form W_S(g) = pole + A_inf - sum_{p in S} prime_p(g), g = f*f~,
# in-tree normalisations (ArithmeticReflectionPositivity / ArchimedeanFloor / WeilSupportLadder).
# f broad bump with int f e^{u/2} = 0 (pole term = 2 A+ A- = 0). Expect W_S/g(0) -> C_inf - 2 sum_S kappa_p.
import numpy as np
from math import log, sqrt, pi
gamma=0.5772156649015329
Cinf=-gamma-log(pi)-pi/2-3*log(2)
def bump(x):
    y=np.zeros_like(x); m=np.abs(x)<1; y[m]=np.exp(-1/(1-x[m]**2)); return y
def run(sig, S, n=8001):
    u=np.linspace(-sig,sig,n); du=u[1]-u[0]
    b=bump(u/sig); w=u/sig*b
    lam=np.sum(b*np.exp(u/2))/np.sum(w*np.exp(u/2)); f=b-lam*w
    Ap=np.sum(f*np.exp(u/2))*du; Am=np.sum(f*np.exp(-u/2))*du
    full=np.correlate(f,f,mode='full')*du
    lags=np.arange(-(n-1),n)*du
    G=full[::-1]
    g=lambda x: np.interp(x,lags,G,left=0.0,right=0.0)
    g0=float(g(0.0))
    uu=np.linspace(1e-7,2*sig+40,400001); gu=g(uu)
    integrand=(np.exp(-2*uu)*g0-np.exp(-uu/2)*gu)/(1-np.exp(-2*uu))
    Ainf=-(gamma+log(pi))*g0+2*np.trapezoid(integrand,uu)
    pole=2*Ap*Am
    prime=0.0
    for p in S:
        k=1
        while k*log(p)<=2*sig:
            prime+=log(p)*p**(-k/2)*(g(k*log(p))+g(-k*log(p))); k+=1
    W=pole+Ainf-prime
    kap=sum(log(p)/(sqrt(p)-1) for p in S)
    return g0,pole/g0,Ainf/g0,W/g0,Cinf-2*kap
if __name__=="__main__":
    print("C_inf =",Cinf)
    for S in ([],[2],[2,3],[2,3,5,7]):
        for sig in (0.3,2,5,10,20,40):
            g0,pl,ai,W,lim=run(sig,S)
            print(f"S={S} sigma={sig}: pole/g0={pl:.1e} Ainf/g0={ai:.4f} W_S/g0={W:.4f} limit={lim:.4f}",flush=True)
