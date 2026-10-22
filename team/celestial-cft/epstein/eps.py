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

def Kbest(z, nu):
    # Arb's K_nu(z) is rigorous, but its radius is erratic near its internal algorithm switch (seen for z ~ 40-140).
    # Every enclosure is valid, so try bessel_k and K = sqrt(pi) (2z)^nu e^{-z} U(nu+1/2, 2nu+1, 2z) at several
    # precisions and keep the tightest one.
    p0=ctx.prec; best=None
    try:
        for p in (p0, p0+64, 2*p0, 4*p0):
            ctx.prec=p; z_=acb(z); nu_=acb(nu)
            for meth in (0,1):
                v = z_.bessel_k(nu_) if meth==0 else acb.pi().sqrt()*(2*z_)**nu_*(-z_).exp()*(2*z_).hypgeom_u(nu_+0.5,2*nu_+1)
                if not v.is_finite(): continue
                r=float(abs(v).rad())
                if best is None or r<best[0]: best=(r,v)
            if best is not None and best[0] < 2.0**(-p0+8)*float(abs(best[1]).mid()): break
    finally:
        ctx.prec=p0
    return +best[1]

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
        S+= acb(n)**(s-half)*sg*Kbest(2*pi*n*y,nu)*(2*pi*n*x).cos()
    # rigorous tail
    mu=arb(max(abs(sig_hi-0.5),abs(sig_lo-0.5)))
    A=sig_hi+max(0.0,1-2*sig_lo)
    M=N+1
    rate=arb(A)/M-2*pi*y
    assert float(rate.mid())<-0.1
    Kb=Kbest(2*pi*M*y,acb(mu)).real
    tail=2*Kb*arb(M)**arb(A)/(1-rate.exp())
    # (tail bounds the bracketed sum S; pref multiplies S below. An earlier version multiplied by |pref| twice.)
    S=S+acb(arb(0,tail.upper()), arb(0,tail.upper()))
    return out+pref*S

def Z_Q(Q, s, N=None):
    # tau_Q must be formed at the raised working precision; otherwise its rounding radius is amplified by ~e^{pi t/2}.
    s=acb(s); t0=abs(float(s.imag.mid())); old=ctx.prec
    ctx.prec=64+int(2.6*t0)+20+old
    try:
        (x,y),D=tau_of(Q)
        pre=(acb(2)/arb(-D).sqrt())**s
        ctx.prec=old
        E=eisenstein(x,y,s,N)
        ctx.prec=64+int(2.6*t0)+20+old
        return _fin(pre*E, old)
    finally:
        ctx.prec=old

def _fin(v, p):
    ctx.prec=p; return +v

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
