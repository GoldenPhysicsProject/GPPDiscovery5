"""C.15 even (K3) certificate framework: float exploration (NUMERICAL) of the Schur bound with the comparison form
Q'_T = q_T H + A_T (Thm C.15.A), monotonicity table, full-Q Galerkin cross-check, and Lemma-A tail-kernel check.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import numpy as np, mpmath as mp, sys
from scipy.special import digamma
a=np.log(2); c2=a/np.sqrt(2)
def Om(t):  # vectorised Re psi(1/4+it/2) - log pi
    return np.real(digamma(0.25+0.5j*np.asarray(t,dtype=complex)))-np.log(np.pi)
Om0=Om(0.0)
def ehat(k,t,L):
    w=(k+0.5)*np.pi/L
    return (L*np.sinc((w-t)*L/np.pi)+L*np.sinc((w+t)*L/np.pi))/np.sqrt(L)
def gl(a0,b0,n):
    x,w=np.polynomial.legendre.leggauss(n); return 0.5*(b0-a0)*(x+1)+a0, 0.5*(b0-a0)*w
def panels(T,width=2.0,n=24):
    ts=[];ws=[]
    edges=np.arange(0,T+1e-12,width)
    if edges[-1]<T: edges=np.append(edges,T)
    for i in range(len(edges)-1):
        t,w=gl(edges[i],edges[i+1],n); ts.append(t); ws.append(w)
    return np.concatenate(ts),np.concatenate(ws)
def Bmat(N,L):
    # B_jk = int_{-L}^{L-a} [e_j'(x) e_k'(x+a) + 1/4 e_j(x) e_k(x+a)] dx, e_k = cos(w_k x)/sqrt(L); Gauss-Legendre exact enough (smooth)
    x,wq=gl(-L,L-a,400); W=(np.arange(N)+0.5)*np.pi/L
    E=np.cos(np.outer(W,x))/np.sqrt(L); Ea=np.cos(np.outer(W,x+a))/np.sqrt(L)
    D=-W[:,None]*np.sin(np.outer(W,x))/np.sqrt(L); Da=-W[:,None]*np.sin(np.outer(W,x+a))/np.sqrt(L)
    B=(D*wq)@Da.T+0.25*(E*wq)@Ea.T
    return 0.5*(B+B.T)
def fullQ(N,L,X=4000.0):
    W=(np.arange(N)+0.5)*np.pi/L; nrm=np.sqrt(W**2+0.25)
    t,wt=panels(X,1.0,16)
    E=np.array([ehat(k,t,L) for k in range(N)])
    R=(E*(wt*(0.25+t*t)*Om(t)/np.pi))@E.T
    # tail beyond X: (1/pi) int_X^inf t^2 Om(t) * 4 w_j w_k (-1)^{j+k} cos^2(tL)/(L t^4) ~ cos^2 -> 1/2
    sg=(-1.0)**np.arange(N)
    tail=(2/(np.pi*L))*(np.log(X/2)-np.log(np.pi)+1)/X   # int_X^inf Om/t^2 ~ (log(X/2pi)+1)/X
    R+=np.outer(sg*W,sg*W)*tail
    Q=R-2*c2*Bmat(N,L)
    return Q/np.outer(nrm,nrm), R/np.outer(nrm,nrm)
def Qprime(N,L,T,width=1.0):
    W=(np.arange(N)+0.5)*np.pi/L; nrm=np.sqrt(W**2+0.25)
    t,wt=panels(T,width,16)
    q=Om(t)-2*c2*np.cos(a*t); qT=Om(T)-2*c2
    E=np.array([ehat(k,t,L) for k in range(N)])
    A=(E*(wt*(0.25+t*t)*(q-qT)/np.pi))@E.T/np.outer(nrm,nrm)
    return A+qT*np.eye(N), qT, Om(T)-Om0
def weighted_tail(N,L,T,wfun,K=6000):
    t,wt=panels(T,1.0,16); g=wt*(0.25+t*t)*wfun(t)/np.pi
    ks=np.arange(N,K); W=(ks+0.5)*np.pi/L
    s=0.0
    for i0 in range(0,len(ks),500):
        kk=ks[i0:i0+500]; WW=W[i0:i0+500]
        E=(L*np.sinc((WW[:,None]-t)*L/np.pi)+L*np.sinc((WW[:,None]+t)*L/np.pi))/np.sqrt(L)
        s+=np.sum((E**2)@g/(WW**2+0.25))
    wK=(K+0.5)*np.pi/L; Wmax=np.max(wfun(t))
    rem=Wmax*(T/4+T**3/3)/np.pi*4/L*(L/np.pi)/(3*(wK-np.pi/L-T)**3)  # rigorous: ehat_k^2/(w^2+1/4) <= 4/(L (w-T)^4), sum <= integral
    return s+rem
def cert(ell,T,N):
    L=(a+ell)/2
    A,qT,DT=Qprime(N,L,T); muA=np.linalg.eigvalsh(A)[0]
    W=(np.arange(N)+0.5)*np.pi/L; nrm=np.sqrt(W**2+0.25)
    t,wt=panels(T,1.0,16); q=Om(t)-2*c2*np.cos(a*t); aa=q-qT
    E=np.array([ehat(k,t,L) for k in range(N)])
    absA=(E*(wt*(0.25+t*t)*np.abs(aa)/np.pi))@E.T/np.outer(nrm,nrm); sV=np.linalg.eigvalsh(absA)[-1]
    absfun=lambda tt: np.abs(Om(tt)-2*c2*np.cos(a*tt)-qT)
    negfun=lambda tt: np.maximum(0,-(Om(tt)-2*c2*np.cos(a*tt)-qT))
    epsAbs=weighted_tail(N,L,T,absfun); epsNeg=weighted_tail(N,L,T,negfun); eps=weighted_tail(N,L,T,lambda tt: np.ones_like(tt))
    B2=sV*epsAbs; ctail=qT-epsNeg
    s=muA+ctail; lam=(s-np.sqrt((muA-ctail)**2+4*B2))/2
    return dict(muA=muA,qT=qT,DT=DT,sV=sV,eps=eps,epsAbs=epsAbs,epsNeg=epsNeg,B=np.sqrt(B2),ctail=ctail,lam=lam)
def KW(t,L,N):   # tail kernel sum_{k>=N} |f_k^(t)|^2 = K_e(t) - sum_{k<N}, K_e from S&M Lemma A (7906287)
    t=np.asarray(t,float); s2=np.where(t>0,np.sin(2*t*L)/(2*np.where(t>0,t,1)),L)
    dens=(L+s2-2*np.cos(t*L)*(0.5*np.tanh(L/2)*np.cos(t*L)+t*np.sin(t*L))/(t*t+0.25))/np.pi
    W=(np.arange(N)+0.5)*np.pi/L
    return np.pi*dens/(0.25+t*t)-sum(ehat(k,t,L)**2/(W[k]**2+0.25) for k in range(N))
if __name__=='__main__':
    print("== Thm C.15.A Schur certificate (float), T=60, N=64 ==")
    for ell in (0.05,0.1,0.15,0.2):
        r=cert(ell,60.0,64)
        print(f"ell={ell}: muA={r['muA']:.5f} q_T={r['qT']:.5f} s_V={r['sV']:.4f} eps^1={r['eps']:.3e} eps^|a|={r['epsAbs']:.3e} eps^-={r['epsNeg']:.3e} ||B||<={r['B']:.5f} c_tail>={r['ctail']:.5f} lam_cert={r['lam']:+.5f}",flush=True)
    print("== N, T dependence at ell=0.1 ==")
    for T,N in ((60.,40),(60.,96),(100.,64)):
        r=cert(0.1,T,N); print(f"T={T} N={N}: muA={r['muA']:.5f} ||B||<={r['B']:.5f} c_tail>={r['ctail']:.4f} lam_cert={r['lam']:+.5f}",flush=True)
    print("== full-Q Galerkin (Dirichlet, N=64), monotone in ell ==")
    for ell in (0.028,0.05,0.1,0.15,0.2):
        Q,R=fullQ(64,(a+ell)/2); print(f"ell={ell}: lam_min(Q)={np.linalg.eigvalsh(Q)[0]:+.5f} lam_min(R)={np.linalg.eigvalsh(R)[0]:+.5f}")
    print("== tail-kernel check: (1/pi) int_0^T (1/4+t^2) K_W dt (Lemma A closed form) vs mode sum eps^1, ell=0.1, T=60, N=64 ==")
    L=(a+0.1)/2; t,wt=panels(60.0,1.0,16)
    print(np.sum(wt*(0.25+t*t)*KW(t,L,64))/np.pi, weighted_tail(64,L,60.0,lambda tt: np.ones_like(tt)))
def cert_majorant(ell,T,N):
    """Same as cert, with the analytic majorants m_abs = Om(T)-Om(t)+2c2(1-cos at) >= |a_T| and m_neg = Om(T)-Om(t) >= (a_T)_- on [0,T]."""
    L=(a+ell)/2
    A,qT,DT=Qprime(N,L,T); muA=np.linalg.eigvalsh(A)[0]
    W=(np.arange(N)+0.5)*np.pi/L; nrm=np.sqrt(W**2+0.25)
    t,wt=panels(T,1.0,16); OmT=Om(T)
    E=np.array([ehat(k,t,L) for k in range(N)])
    mabs=lambda tt: OmT-Om(tt)+2*c2*(1-np.cos(a*tt)); mneg=lambda tt: OmT-Om(tt)
    sV=np.linalg.eigvalsh((E*(wt*(0.25+t*t)*mabs(t)/np.pi))@E.T/np.outer(nrm,nrm))[-1]
    eA=weighted_tail(N,L,T,mabs); eN=weighted_tail(N,L,T,mneg)
    B2=sV*eA; ct=qT-eN; lam=(muA+ct-np.sqrt((muA-ct)**2+4*B2))/2
    return dict(muA=muA,sV=sV,epsAbs=eA,epsNeg=eN,B=np.sqrt(B2),ctail=ct,lam=lam)
if __name__=='__main__':
    print("== analytic-majorant version (recommended for Arb: all integrands analytic on [0,T]) ==")
    for ell in (0.1,0.2):
        r=cert_majorant(ell,60.0,64)
        print(f"ell={ell} T=60 N=64: muA={r['muA']:.5f} s_V^m={r['sV']:.4f} eps^m_abs={r['epsAbs']:.3e} eps^m_neg={r['epsNeg']:.3e} ||B||<={r['B']:.5f} c_tail>={r['ctail']:.5f} lam_cert={r['lam']:+.5f}",flush=True)
