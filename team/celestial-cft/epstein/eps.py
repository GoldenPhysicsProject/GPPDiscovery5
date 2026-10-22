"""Rigorous (ball-arithmetic) evaluation of real-analytic Eisenstein series / Epstein zeta.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

E(tau,s) = sum'_{(m,n)} y^s/|m tau+n|^{2s}
 = 2 zeta(2s) y^s + 2 sqrt(pi) Gamma(s-1/2)/Gamma(s) zeta(2s-1) y^{1-s}
   + 8 pi^s sqrt(y)/Gamma(s) sum_{n>=1} n^{s-1/2} sigma_{1-2s}(n) K_{s-1/2}(2 pi n y) cos(2 pi n x)
Tail beyond N is enclosed rigorously using |K_nu(z)| <= K_{Re nu}(z) (z>0) and the
monotone decrease of e^z K_mu(z); |sigma_{1-2s}(n)| <= d(n) n^{max(0,1-2 sigma)} <= 2 sqrt(n) n^{max(0,1-2sigma)}.
For a form Q=[a,b,c], disc D=b^2-4ac<0, tau_Q=(-b+i sqrt|D|)/(2a):
  E(tau_Q,s) = (sqrt|D|/2)^s Z_Q(s),  Z_Q(s)=sum' Q(m,n)^{-s};
  Lambda_Q(s) = (sqrt|D|/(2pi))^s Gamma(s) Z_Q(s) = pi^{-s} Gamma(s) E(tau_Q,s) = Lambda_Q(1-s).
"""
from flint import acb, arb, ctx
import math

def divisors(n):
    r=[]; d=1
    while d*d<=n:
        if n%d==0:
            r.append(d)
            if d*d!=n: r.append(n//d)
        d+=1
    return r

def tau_of(Q):
    a,b,c=Q; D=b*b-4*a*c
    return (arb(-b)/(2*a), arb(-D).sqrt()/(2*a)), D

def eisenstein(x, y, s, N=None, base=64):
    s=acb(s)
    t0=abs(float(s.imag.mid()))
    old=ctx.prec
    ctx.prec=base+int(2.6*t0)+20
    try:
        return _eis(x,y,s,N,old)
    finally:
        ctx.prec=old

def _eis(x,y,s,N,outprec):
    pi=arb.pi()
    sig_hi = float((s.real.mid()+s.real.rad()))
    sig_lo = float((s.real.mid()-s.real.rad()))
    t=abs(float(s.imag.mid()))+float(s.imag.rad())
    yf=float(y.mid())
    if N is None:
        N=int((math.pi*t/2 + 0.75*outprec + 10)/(2*math.pi*yf))+2
    half=acb(0.5)
    out = 2*(2*s).zeta()*acb(y)**s + 2*pi.sqrt()*(s-half).gamma()*s.rgamma()*(2*s-1).zeta()*acb(y)**(1-s)
    pref = 8*acb(pi)**s*y.sqrt()*s.rgamma()
    nu=s-half
    S=acb(0)
    for n in range(1,N+1):
        sg=acb(0)
        for d in divisors(n): sg+= acb(d)**(1-2*s)
        S+= acb(n)**(s-half)*sg*acb(2*pi*n*y).bessel_k(nu)*(2*pi*n*x).cos()
    # rigorous tail
    mu=arb(max(abs(sig_hi-0.5),abs(sig_lo-0.5)))
    A=sig_hi+max(0.0,1-2*sig_lo)
    M=N+1
    rate=arb(A)/M-2*pi*y
    assert float(rate.mid())<-0.1
    Kb=acb(2*pi*M*y).bessel_k(acb(mu)).real
    tail=2*Kb*arb(M)**arb(A)/(1-rate.exp())
    tail=abs(pref)*tail
    S=S+acb(arb(0,tail.upper()), arb(0,tail.upper()))
    return out+pref*S

def Z_Q(Q, s, N=None):
    (x,y),D=tau_of(Q)
    s=acb(s)
    return (acb(2)/arb(-D).sqrt())**s*eisenstein(x,y,s,N)

def Lambda_Q(Q, s, N=None):
    (x,y),D=tau_of(Q)
    s=acb(s)
    return acb(arb.pi())**(-s)*s.gamma()*eisenstein(x,y,s,N)

def Z_direct(Q,s,M=400):
    a,b,c=Q; tot=0
    for m in range(-M,M+1):
        for n in range(-M,M+1):
            if m==0 and n==0: continue
            tot+= (a*m*m+b*m*n+c*n*n)**(-s)
    return tot
