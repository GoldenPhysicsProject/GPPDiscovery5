"""NUMERICAL (float64, independent of s27): locate lambda_min and the extremal test of the zeta Weil form
on supp f in [-h, h] (in-tree rung window L_tree = 2h).  No zeros used.
(a) Galerkin of Zhu's lower form R_T (Q >= R_T, R_T increases to Q as T -> inf) on orthonormal Legendre modes:
    R_T(f) = 2<cosh(x/2),f>^2 - 2<sinh(x/2),f>^2 + (1/pi) int_0^T (Psi-beta)|F|^2 dt + beta||f||^2,
    Psi(t) = Re psi(1/4+it/2) - log pi - sum_{log n<2h} 2 Lambda(n)/sqrt(n) cos(t log n), beta = log(T/2pi) - 1/T - A.
(b) Rayleigh quotient of the TRUE form Q at the extremal vector, in the time domain with the in-tree
    definitions (poleTerm, archimedeanTerm via the Levy split A = C_inf g(0) + 2 int e^{-u/2}(g(0)-g(u))/(1-e^{-2u}),
    primeSide), giving an upper bound for lambda_min(Q)."""
import numpy as np, math, sys, json
from scipy.special import spherical_jn, psi, eval_legendre
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss
gamma=0.5772156649015329; Cinf=-gamma-math.log(math.pi)-math.pi/2-3*math.log(2)
def mangoldt(n):
    for p in range(2,n+1):
        if n%p==0:
            m=n
            while m%p==0: m//=p
            return (math.log(p) if m==1 else 0.0)
def comb(h):
    return [(n,mangoldt(n)) for n in range(2,int(math.exp(2*h))+2) if math.log(n)<2*h and mangoldt(n)>0]
def run(h,K=40,T=800.0,nt=200000):
    terms=comb(h); A=sum(2*l/math.sqrt(n) for n,l in terms)
    beta=math.log(T/(2*math.pi))-1/T-A
    xg,wg=leggauss(400); x=h*xg; w=h*wg
    nrm=np.array([math.sqrt((2*k+1)/(2*h)) for k in range(K)])
    P=np.array([eval_legendre(k,xg) for k in range(K)])*nrm[:,None]   # orthonormal on [-h,h]
    c=P@(w*np.cosh(x/2)); s=P@(w*np.sinh(x/2))
    g20,w20=leggauss(20); npan=int(T/0.5); a=np.arange(npan)*0.5
    t=(a[:,None]+(g20[None,:]+1)*0.25).ravel(); tw=np.tile(w20*0.25,npan)
    Psi=psi(0.25+0.5j*t).real-math.log(math.pi)-sum(2*l/math.sqrt(n)*np.cos(t*math.log(n)) for n,l in terms)
    J=np.array([spherical_jn(k,h*t) for k in range(K)])*nrm[:,None]*2*h   # |F_k(t)| up to phase i^k
    out={}
    for par in (0,1):
        idx=[k for k in range(K) if k%2==par]
        # Re(conj F_j F_k) = Re(i^{k-j}) J_j J_k ; same parity => i^{k-j} = (-1)^{(k-j)/2}
        sg=np.array([[(-1)**((k-j)//2) for k in idx] for j in idx])
        M=(J[idx]*(tw*(Psi-beta))) @ J[idx].T * sg / math.pi
        R=2*np.outer(c[idx],c[idx])-2*np.outer(s[idx],s[idx])+M+beta*np.eye(len(idx))
        ev,V=np.linalg.eigh(R); out[par]=(ev[0],V[:,0],idx)
    return out,P,x,w,terms,nrm
def Qtrue(h,coef,idx,terms,K):
    nrm=np.array([math.sqrt((2*k+1)/(2*h)) for k in range(K)])
    def f(y):
        y=np.atleast_1d(y); r=np.zeros_like(y,dtype=float); m=np.abs(y)<=h
        r[m]=sum(cf*nrm[k]*eval_legendre(k,y[m]/h) for cf,k in zip(coef,idx)); return r
    xg,wg=leggauss(300)
    def g(u):
        u=abs(u)
        if u>=2*h: return 0.0
        a,b=-h+u,h; y=(xg+1)*(b-a)/2+a; return float(np.sum(wg*(b-a)/2*f(y)*f(y-u)))
    g0=g(0.0)
    Ap=float(np.sum(wg*h*f(h*xg)*np.exp(h*xg/2))); Am=float(np.sum(wg*h*f(h*xg)*np.exp(-h*xg/2)))
    pole=2*Ap*Am
    I1,e1=quad(lambda u: math.exp(-u/2)*(g0-g(u))/(1-math.exp(-2*u)),0,2*h,limit=400,epsabs=1e-15,epsrel=1e-13)
    I2=g0*quad(lambda u: math.exp(-u/2)/(1-math.exp(-2*u)),2*h,80,epsabs=1e-16,epsrel=1e-14)[0]
    arch=Cinf*g0+2*(I1+I2)
    prime=sum(l/math.sqrt(n)*2*g(math.log(n)) for n,l in terms)
    return (pole+arch-prime)/g0, dict(pole=pole/g0,arch=arch/g0,prime=prime/g0,g_log2=(g(math.log(2))/g0))
if __name__=="__main__":
    res={}
    for name,h in (("rung0_edge",0.3465),("rung1_mid",0.45),("rung1_edge",0.5493),("rung2_edge",0.6931)):
        for T in (200.0,800.0,3200.0):
            out,P,x,w,terms,nrm=run(h,K=40,T=T)
            print(name,"h=",h,"T=",T,"lam_min(R_T) even/odd:",f"{out[0][0]:.4e} {out[1][0]:.4e}",flush=True)
        par=0 if out[0][0]<out[1][0] else 1
        lam,v,idx=out[par]
        q,parts=Qtrue(h,v,idx,terms,40)
        prof=[float(sum(cf*nrm[k]*eval_legendre(k,xx/h) for cf,k in zip(v,idx))) for xx in np.linspace(-h,h,11)]
        print(name,"extremal parity",["even","odd"][par],"Q(v)/|v|^2 =",f"{q:.6e}",parts,flush=True)
        print(name,"profile on 11 pts:",np.round(np.array(prof)/max(abs(np.array(prof))),4).tolist(),flush=True)
        res[name]=dict(h=h,comb=[n for n,_ in terms],lamR_T3200=[float(out[0][0]),float(out[1][0])],Q_at_extremal=q,parts=parts,parity=["even","odd"][par],profile=prof)
    json.dump(res,open("rung_extremal.json","w"),indent=1)
