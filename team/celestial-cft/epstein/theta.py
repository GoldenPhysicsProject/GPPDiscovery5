"""Rigorous Lambda_Q(s) via theta-function splitting (no Bessel functions).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

For positive definite Q=[a,b,c], D=b^2-4ac<0, v(m,n)=2 pi Q(m,n)/sqrt|D|:
  Lambda_Q(s) := (sqrt|D|/(2 pi))^s Gamma(s) Z_Q(s)
             = 1/(s-1) - 1/s + sum'_{(m,n)} [G(s,v)+G(1-s,v)],  G(s,v)=v^{-s} Gamma(s,v).
(Poisson: theta_Q(1/x)=x theta_{Q*}(x), Q*=[c,-b,a] has the same values as Q.)
Tail over Q(m,n)>Qmax enclosed with |G(s,v)| <= e^{-v}/(v-a), a=max(0,Re s-1) and the
lattice count #{Q<=X} <= pi (sqrt(X/lam)+sqrt 2)^2, lam = min eigenvalue of [[a,b/2],[b/2,c]].
"""
from flint import acb, arb, ctx
import math
from functools import lru_cache

@lru_cache(None)
def values(Q, Qmax):
    a,b,c=Q
    lam=(a+c-math.sqrt((a-c)**2+b*b))/2
    R=int(math.sqrt(Qmax/lam))+2
    cnt={}
    for m in range(-R,R+1):
        for n in range(-R,R+1):
            if m==0 and n==0: continue
            q=a*m*m+b*m*n+c*n*n
            if q<=Qmax: cnt[q]=cnt.get(q,0)+1
    return tuple(sorted(cnt.items())), lam

def tail_bound(Q, Qmax, sig_hi, sig_lo):
    a,b,c=Q; D=b*b-4*a*c
    vals,lam=values(Q,Qmax)
    cst=2*math.pi/math.sqrt(-D)
    A=max(0.0, sig_hi-1, -sig_lo)  # covers both G(s,.) and G(1-s,.)
    tot=arb(0); k=Qmax
    while True:
        Nk=arb.pi()*(arb(k+1)/lam).sqrt().__add__(arb(2).sqrt())**2
        v=arb(cst)*k
        term=2*Nk*(-v).exp()/(v-A)
        tot+=term
        if float(term.mid())<1e-60*max(float(tot.mid()),1e-300) and k>Qmax+50:
            # geometric remainder (ratio <= e^{-cst}*(1+...)<=e^{-cst/2})
            r=math.exp(-cst/2)
            tot+=term*arb(r)/(1-arb(r))
            break
        k+=1
    return tot

def Lambda(Q, s, bits=60):
    """bits: target relative accuracy; picks Qmax from height."""
    s=acb(s)
    a,b,c=Q; D=b*b-4*a*c
    t=abs(float(s.imag.mid()))+float(s.imag.rad())
    sig_hi=float(s.real.mid()+s.real.rad()); sig_lo=float(s.real.mid()-s.real.rad())
    vmax=math.pi*t/2+bits*math.log(2)+20+2*math.log(t+3)
    Qmax=int(vmax*math.sqrt(-D)/(2*math.pi))+1
    vals,_=values(Q,Qmax)
    old=ctx.prec; ctx.prec=max(old, bits+int(2.3*t)+40)
    try:
        cst=2*arb.pi()/arb(-D).sqrt()
        one=acb(1)
        out=one/(s-1)-one/s
        sp=one-s
        for q,mlt in vals:
            v=acb(cst*q)
            out+= mlt*(v**(-s)*v.gamma_upper(s) + v**(-sp)*v.gamma_upper(sp))
        tb=tail_bound(Q,Qmax,sig_hi,sig_lo)
        out+= acb(arb(0,tb.upper()),arb(0,tb.upper()))
        return +out
    finally:
        ctx.prec=old
