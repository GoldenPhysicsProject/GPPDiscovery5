# NUMERICAL (float) probe for P.17: full {2,3} Weil-form minima vs 2L, and a C.15.A-type Schur bound at the right endpoint.
import numpy as np
from scipy.special import psi
from scipy.linalg import eigh
a=np.log(2); b=np.log(3); c2=a/np.sqrt(2); c3=b/np.sqrt(3)
Om=lambda t: psi(0.25+0.5j*t).real-np.log(np.pi)
def ft(n,t,L):
    k=n*np.pi/(2*L)
    I=lambda w: 2*L*np.sinc(w*L/np.pi)
    return (np.exp(1j*k*L)*I(k-t)-np.exp(-1j*k*L)*I(-k-t))/(2j)
def dens(par,t,L):
    if par=='odd':
        return (L - L*np.sinc(2*t*L/np.pi) - 2*np.sin(t*L)*(0.5/np.tanh(L/2)*np.sin(t*L)-t*np.cos(t*L))/(t*t+.25))/np.pi
    return (L + np.sin(2*t*L)/(2*t) - 2*np.cos(t*L)*(0.5*np.tanh(L/2)*np.cos(t*L)+t*np.sin(t*L))/(t*t+.25))/np.pi
def modes(par,N,L,t):
    ns=[n for n in range(1,4*N+4) if (n%2==0)==(par=='odd')][:N]
    F=np.array([ft(n,t,L)/np.sqrt(((n*np.pi/(2*L))**2+.25)*L) for n in ns])  # H-normalised
    return F
def full_min(par,L,N=60,Tmax=3000.,ng=300000):
    t=np.linspace(0,Tmax,ng+1)[1:]; dt=t[1]-t[0]; F=modes(par,N,L,t); w=(0.25+t*t)/np.pi*dt
    s=Om(t)-2*c2*np.cos(a*t)-(2*c3*np.cos(b*t) if 2*L>b else 0)
    M=np.real((F*(w*s))@F.conj().T); G=np.real((F*w)@F.conj().T)
    return eigh(M,G,eigvals_only=True)[0], np.abs(G-np.eye(N)).max()
def schur(par,L,T,N,ng=200000):
    t=np.linspace(1e-7,T,ng+1); dt=t[1]-t[0]; F=modes(par,N,L,t); w=(0.25+t*t)/np.pi*dt
    w[0]*=.5; w[-1]*=.5
    lines=lambda x: 2*c2*(1-np.cos(a*x))+(2*c3*(1-np.cos(b*x)) if 2*L>b else 0)
    qT=Om(T)-2*c2-(2*c3 if 2*L>b else 0)
    aT=Om(t)-Om(T)+lines(t)
    mabs=Om(T)-Om(t)+lines(t); mneg=Om(T)-Om(t)
    gram=lambda g: np.real((F*(w*g))@F.conj().T)
    A=qT*np.eye(N)+gram(aT); muA=eigh(A,eigvals_only=True)[0]
    sV=eigh(gram(mabs),eigvals_only=True)[-1]
    KW=dens(par,t,L)*np.pi/(0.25+t*t)-np.sum(np.abs(F)**2,axis=0)
    eps=lambda m: np.sum(w*m*KW)
    ea,en=eps(mabs),eps(mneg); beta=np.sqrt(max(sV*ea,0)); ct=qT-en
    lam=0.5*(muA+ct-np.sqrt((muA-ct)**2+4*beta**2))
    return dict(qT=qT,muA=muA,sV=sV,eabs=ea,eneg=en,beta=beta,ctail=ct,lam=lam)
if __name__=='__old__':
    for par in ['odd','even']:
        for twoL in [b, np.log(3.5), np.log(3.9), np.log(4)-1e-9]:
            m,go=full_min(par,twoL/2)
            print('%s 2L=%.5f full{2,3} Ritz min <= %.3e (Gram dev %.1e)'%(par,twoL,m,go))

def main_runs():
    print('# p17_probe.py -- NUMERICAL (float64, trapezoid quadrature, no remainders). Not a certificate. RH OPEN.')
    print('# Part A: Ritz upper bounds for the full {2,3} Weil form (sine basis N=60, Tmax=3000).')
    for par in ['odd','even']:
        for l3 in [0.0,0.01,0.03,0.05,0.1]:
            m,go=full_min(par,(b+l3)/2)
            print('A %s l3=%.2f min<=%.4e gramdev=%.1e'%(par,l3,m,go),flush=True)
        for twoL in [np.log(3.5), np.log(3.9), np.log(4)-1e-9]:
            m,go=full_min(par,twoL/2)
            print('A %s 2L=%.5f min<=%.4e gramdev=%.1e'%(par,twoL,m,go),flush=True)
    print('# Part B: C.15.A-type Schur bound extended to {2,3} (float, Riemann-sum eps).')
    for par,l3,T,N in [('odd',0.01,100,64),('odd',0.01,100,128),('odd',0.01,100,200),('odd',0.01,60,128),
                       ('odd',0.01,100,300),('odd',0.02,100,300),('odd',0.03,100,300),('odd',0.05,100,300),('even',0.0,100,300)]:
        r=schur(par,(b+l3)/2,T,N,ng=40000)
        print('B',par,l3,T,N,' '.join('%s=%.4e'%(k,v) for k,v in r.items()),flush=True)

if __name__=='__main__':
    main_runs()
