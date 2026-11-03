# NUMERICAL probe (not a certificate): is (H3^o) at 2L=log3, T=16, eta=0.0373 even plausible?
import numpy as np
from scipy.special import psi
from scipy.linalg import eigh
a=np.log(2); b=np.log(3); c2=a/np.sqrt(2); L=b/2
def Om(t): return psi(0.25+0.5j*t).real-np.log(np.pi)
def basis_ft(n,t):  # FT of sin(n pi (x+L)/2L) on (-L,L)
    k=n*np.pi/(2*L)
    # int_{-L}^{L} sin(k(x+L)) e^{-itx} dx
    def I(w): # int_{-L}^{L} e^{i w x} dx
        return np.where(np.abs(w)<1e-12, 2*L, 2*np.sin(w*L)/np.where(np.abs(w)<1e-12,1,w))
    return (np.exp(1j*k*L)*I(k-t)-np.exp(-1j*k*L)*I(-k-t))/(2j)
def run(par,N,T,Tmax=4000.0,ng=400000):
    ns=[n for n in range(1,4*N) if (n%2==0)==(par=='odd')][:N]
    t=np.linspace(0,Tmax,ng+1)[1:]; dt=t[1]-t[0]
    F=np.array([basis_ft(n,t) for n in ns])
    w=(0.25+t*t)/np.pi*dt
    def gram(sym): return np.real((F*(w*sym))@F.conj().T)
    G=gram(np.ones_like(t)); hi=(t>=T).astype(float)
    Q2=gram(Om(t)-2*c2*np.cos(a*t)); Nhi=gram(hi); Nlo=G-Nhi
    m2=eigh(Q2,G,eigvals_only=True)[0]
    H0=Q2-Om(T)*Nhi
    lam=eigh(H0,G,eigvals_only=True)[0]
    # best mu with H0 - mu*Nlo >= 0 : bisection
    lo,hi_=-1.0,1.0
    for _ in range(60):
        mid=(lo+hi_)/2
        if eigh(H0-mid*Nlo,G,eigvals_only=True)[0]>=0: lo=mid
        else: hi_=mid
    return m2,lam,lo
for par in ['odd','even']:
    for N in [20,40,60]:
        m2,lam,mu=run(par,N,16.0)
        print(par,N,'P2 margin<=%.5f'%m2,' lambda_min(Q2-Om16*Nhi)<=%.5f'%lam,' mu_max(H3,T=16,eta=.0373)<=%.5f'%mu)
