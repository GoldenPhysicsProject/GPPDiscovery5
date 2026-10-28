"""Grok audit S103i(2) of Celestial C.15 (NUMERICAL cross-checks; run from repo root). Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import numpy as np, mpmath as mp, sys
sys.path.insert(0,'team/celestial-cft/arp_c')
import c15_even_schur as S
mp.mp.dps=25
a=mp.log(2); ell=mp.mpf('0.1'); L=(a+ell)/2
om=lambda k:(k+mp.mpf(1)/2)*mp.pi/L
# independent B_jk via direct quad
def Bq(j,k):
    f=lambda x: (om(j)*om(k)*mp.sin(om(j)*x)*mp.sin(om(k)*(x+a)) + mp.cos(om(j)*x)*mp.cos(om(k)*(x+a))/4)/L
    return mp.quad(f,[-L,(-L+L-a)/2,L-a])
print("B00,B13,B33 quad:",[mp.nstr(Bq(*p),14) for p in ((0,0),(1,3),(3,3))])
# K_e vs mode sum
Lf=float(L)
for t in (0.0,5.0,37.0):
    ke=float(S.KW(np.array([t]),Lf,0)[0]); ms=sum(S.ehat(k,t,Lf)**2/(((k+.5)*np.pi/Lf)**2+.25) for k in range(200000))
    print("K_e t=%g closed %.12e modesum %.12e"%(t,ke,ms))
# naive split: normalized full R entries vs 1/(2(j+k+1))
Q,R=S.fullQ(48,Lf)
for (j,k) in ((0,1),(3,10),(10,30),(20,40),(5,40)):
    print("Rt[%d,%d]=%+.5f  1/(2(j+k+1))=%.5f"%(j,k,R[j,k],1/(2*(j+k+1))))
# off-diagonal block norm of full R between V_N and complement (within 48)
for N in (8,16,24):
    print("N=%d ||R[:N,N:48]||=%.4f"%(N,np.linalg.norm(R[:N,N:],2)))
# Omega monotone check
t=np.linspace(0,200,200001); print("Omega increasing on grid:",np.all(np.diff(S.Om(t))>0))
# independent A_N, lam_cert at ell=0.1 via a different quadrature (finer panels, 32-pt)
def QpN(N,T=60.0):
    W=(np.arange(N)+.5)*np.pi/Lf; nrm=np.sqrt(W**2+.25)
    t,wt=S.panels(T,0.5,32); q=S.Om(t)-2*S.c2*np.cos(S.a*t); qT=S.Om(T)-2*S.c2
    E=np.array([S.ehat(k,t,Lf) for k in range(N)])
    A=(E*(wt*(.25+t*t)*(q-qT)/np.pi))@E.T/np.outer(nrm,nrm)+qT*np.eye(N); return A
A=QpN(64); print("muA indep:",np.linalg.eigvalsh(A)[0])
r=S.cert_majorant(0.1,60.0,64); print("cert_majorant 0.1:",r)
# epsilon via closed-form K_W in mpmath (high precision, avoids cancellation)
def KWmp(t,N):
    t=mp.mpf(t)
    s2=mp.sin(2*t*L)/(2*t) if t!=0 else L
    Ke=(L+s2-2*mp.cos(t*L)*(mp.tanh(L/2)/2*mp.cos(t*L)+t*mp.sin(t*L))/(t*t+mp.mpf(1)/4))/(mp.mpf(1)/4+t*t)
    sm=0
    for k in range(N):
        w=om(k); e=2*w*mp.cos(t*L)/(mp.sqrt(L)*(w*w-t*t)); sm+=e*e/(w*w+mp.mpf(1)/4)
    return Ke-sm
OmT=mp.re(mp.digamma(mp.mpf(1)/4+30j))-mp.log(mp.pi); c2=a/mp.sqrt(2)
Omm=lambda t: mp.re(mp.digamma(mp.mpf(1)/4+0.5j*t))-mp.log(mp.pi)
mneg=lambda t:(OmT-Omm(t)); mabs=lambda t: OmT-Omm(t)+2*c2*(1-mp.cos(a*t))
mp.mp.dps=30
eN=mp.quad(lambda t:(mp.mpf(1)/4+t*t)*mneg(t)*KWmp(t,64),mp.linspace(0,60,61))/mp.pi
eA=mp.quad(lambda t:(mp.mpf(1)/4+t*t)*mabs(t)*KWmp(t,64),mp.linspace(0,60,61))/mp.pi
print("eps_neg mp:",mp.nstr(eN,8)," eps_abs mp:",mp.nstr(eA,8))
