"""C.16/C.17: general-L, general-prime-set, both-parity version of the Thm C.15.A Schur certificate (float, NUMERICAL).
Q_L(xi) = (1/2pi) int (1/4+t^2) q(t) |xi^(t)|^2 dt,  q(t) = Omega(t) - sum_{n in S} 2 c_n cos(t log n),  c_n = Lambda(n)/sqrt(n),
for real xi in H^1_0(-L,L) (either parity); S = prime powers with log n < 2L (n with log n = 2L contributes exactly 0).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import numpy as np, sys
from scipy.special import digamma
from sympy import factorint

def Om(t): return np.real(digamma(0.25+0.5j*np.asarray(t,dtype=complex)))-np.log(np.pi)
def vonmangoldt(n):
    f=factorint(n); return np.log(list(f)[0]) if len(f)==1 else 0.0
def prime_set(L, include_endpoint=False, tol=1e-12):
    """[(n, c_n)] for prime powers n with log n < 2L (or <= 2L if include_endpoint)."""
    out=[]; n=2
    while np.log(n) < 2*L+tol:
        lam=vonmangoldt(n)
        if lam>0 and (np.log(n) < 2*L-tol or include_endpoint): out.append((n, lam/np.sqrt(n)))
        n+=1
    return out
def omegas(N,L,parity):
    k=np.arange(N)
    return (k+0.5)*np.pi/L if parity=='even' else (k+1)*np.pi/L
def ehat_w(W,t,L,parity):
    """real transform of e_k (L^{-1/2} cos(w x) even / L^{-1/2} sin(w x) odd); |xi^|^2 = ehat^2. Shape (len W, len t)."""
    W=np.asarray(W)[:,None]; sgn=1.0 if parity=='even' else -1.0
    return np.sqrt(L)*(np.sinc((W-t)*L/np.pi)+sgn*np.sinc((W+t)*L/np.pi))
def gl(a0,b0,n):
    x,w=np.polynomial.legendre.leggauss(n); return 0.5*(b0-a0)*(x+1)+a0, 0.5*(b0-a0)*w
def panels(T,width=1.0,n=16):
    e=np.arange(0,T+1e-12,width)
    if e[-1]<T: e=np.append(e,T)
    ts,ws=zip(*[gl(e[i],e[i+1],n) for i in range(len(e)-1)]); return np.concatenate(ts),np.concatenate(ws)
def qsym(t,P): return Om(t)-sum(2*c*np.cos(np.log(n)*t) for n,c in P)
def qT_of(T,P): return Om(T)-2*sum(c for n,c in P)
def mabs(t,T,P): return Om(T)-Om(t)+sum(2*c*(1-np.cos(np.log(n)*t)) for n,c in P)
def mneg(t,T,P): return Om(T)-Om(t)
def band_matrix(N,L,T,parity,g,width=1.0):
    W=omegas(N,L,parity); nrm=np.sqrt(W**2+0.25); t,wt=panels(T,width)
    E=ehat_w(W,t,L,parity)
    return (E*(wt*(0.25+t*t)*g(t)/np.pi))@E.T/np.outer(nrm,nrm)
def weighted_tail(N,L,T,parity,g,K=None):
    """eps^g_N = (1/pi) int_0^T (1/4+t^2) g K_W dt, mode sum k=N..K-1 + rigorous remainder (f_k^2 <= 4/(L(w_k-T)^4), both parities)."""
    t,wt=panels(T); gg=wt*(0.25+t*t)*g(t)/np.pi
    if K is None: K=max(6000,int(40*T*L)+N)
    W=omegas(K,L,parity)[N:]; s=0.0
    for i in range(0,len(W),500):
        WW=W[i:i+500]; s+=np.sum((ehat_w(WW,t,L,parity)**2)@gg/(WW**2+0.25))
    wK=omegas(K+1,L,parity)[-1]
    rem=np.max(g(t))*(T/4+T**3/3)/np.pi*4/L*(L/np.pi)/(3*(wK-np.pi/L-T)**3)
    return s+rem
def Kclosed(t,L,parity):
    """K_sigma(t) = sum_k f_k^(t)^2 over all modes (even: S&M Lemma A; odd: Lemma C.16.K)."""
    t=np.asarray(t,float); tt=np.where(t>0,t,1.0); s2=np.where(t>0,np.sin(2*t*L)/(2*tt),L)
    c,s=np.cos(t*L),np.sin(t*L)
    if parity=='even': num=L+s2-2*c*(0.5*np.tanh(L/2)*c+t*s)/(t*t+0.25)
    else: num=L-s2-2*s*(0.5/np.tanh(L/2)*s-t*c)/(t*t+0.25)
    return num/(0.25+t*t)
def cert(L,P,T,N,parity,majorant='analytic',width=1.0):
    qT=qT_of(T,P)
    A=band_matrix(N,L,T,parity,lambda t: qsym(t,P)-qT,width)+qT*np.eye(N)
    muA=np.linalg.eigvalsh(A)[0]
    if majorant=='analytic': ga=lambda t: mabs(t,T,P); gn=lambda t: mneg(t,T,P)
    else: ga=lambda t: np.abs(qsym(t,P)-qT); gn=lambda t: np.maximum(0,-(qsym(t,P)-qT))
    M=band_matrix(N,L,T,parity,ga,width); sV=np.linalg.eigvalsh(M)[-1]; sF=np.linalg.norm(M)
    eA=weighted_tail(N,L,T,parity,ga); eN=weighted_tail(N,L,T,parity,gn)
    B2=sV*eA; ct=qT-eN; lam=(muA+ct-np.sqrt((muA-ct)**2+4*B2))/2
    lamF=(muA+ct-np.sqrt((muA-ct)**2+4*sF*eA))/2
    # float rehearsal of the Arb steps: Cholesky of A - mu I at mu = 0.999 muA (when muA>0), and of s I - M at s = 1.001 sV
    chol_ok=None
    if muA>0:
        try: np.linalg.cholesky(A-0.999*muA*np.eye(N)); np.linalg.cholesky(1.001*sV*np.eye(N)-M); chol_ok=True
        except np.linalg.LinAlgError: chol_ok=False
    return dict(L=L,T=T,N=N,parity=parity,qT=qT,muA=muA,sV=sV,sF=sF,eA=eA,eN=eN,B=np.sqrt(B2),ctail=ct,lam=lam,lamF=lamF,chol=chol_ok,A=A)
def galerkin_full(L,P,N,parity,X=6000.0):
    """full-Q Galerkin (upper bound for the sector minimum), Dirichlet modes, |t|<X by quadrature + asymptotic tail."""
    W=omegas(N,L,parity); nrm=np.sqrt(W**2+0.25); t,wt=panels(X,1.0,16)
    E=ehat_w(W,t,L,parity); Qm=(E*(wt*(0.25+t*t)*qsym(t,P)/np.pi))@E.T
    # tail |t|>X: ehat_j ehat_k ~ 4 w_j w_k (+-) trig^2(tL)/(L t^4), trig^2 -> 1/2; int_X^inf Om/t^2 ~ (log(X/2)-log pi + 1)/X
    sg=(-1.0)**np.arange(N); Qm+=np.outer(sg*W,sg*W)*(2/(np.pi*L))*(np.log(X/2)-np.log(np.pi)+1)/X
    Qm/=np.outer(nrm,nrm); return np.linalg.eigvalsh(Qm)[0]
def fmt(r): return (f"2L={2*r['L']:.6f} {r['parity']:4s} T={r['T']:.0f} N={r['N']}: q_T={r['qT']:.5f} mu_A={r['muA']:+.6f} s_V={r['sV']:.4f} "
                    f"(Frob {r['sF']:.3f}) eps^abs={r['eA']:.3e} eps^neg={r['eN']:.3e} ||B||<={r['B']:.5f} c_tail>={r['ctail']:.5f} "
                    f"lam_cert={r['lam']:+.6f} (Frob: {r['lamF']:+.6f}) chol={r['chol']}")

if __name__=='__main__':
    a=np.log(2)
    print("== C.16 odd p=2 Schur certificate (Thm C.16.A, float, NUMERICAL); also general-code regression vs C.15 even ==")
    for par in ('even','odd'):
        for t in (0.0,3.0,10.0,40.0):
            W=omegas(20000,0.4,par); ms=np.sum(ehat_w(W,np.array([t]),0.4,par)[:,0]**2/(W**2+.25))
            print(f"kernel check {par} L=0.4 t={t}: closed {Kclosed(np.array([t]),0.4,par)[0]:.15e} mode-sum(20000) {ms:.15e}")
    L=(a+0.1)/2; print("regression even ell=0.1:",fmt(cert(L,prime_set(L),60.,64,'even')))
    for ell in (0.1,0.2,0.3,0.35,0.4):
        L=(a+ell)/2; P=prime_set(L)
        for maj in ('analytic','sharp'):
            print(f"ell={ell} {maj:8s}",fmt(cert(L,P,60.,64,'odd',maj)),flush=True)
        print(f"ell={ell} full-Q odd Galerkin N=64 (upper bound): {galerkin_full(L,P,64,'odd'):+.6f}",flush=True)
    for T,N in ((60.,40),(60.,96),(100.,64),(100.,96)):
        for ell in (0.1,0.2):
            L=(a+ell)/2; print(f"ell={ell} vary T,N:",fmt(cert(L,prime_set(L),T,N,'odd')),flush=True)
