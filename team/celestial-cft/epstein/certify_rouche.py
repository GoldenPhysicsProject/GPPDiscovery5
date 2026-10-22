"""CERTIFIED isolation of the numerically located zeros of Z_Q for D = -23 (class group C3) and D = -56 (C4), where no
genus-character / Dirichlet-L decomposition exists. Method: Rouche on a small disk, with
  (i)  rigorous point values of Z_Q from eps.py (Chowla-Selberg / K-Bessel expansion in Arb with a rigorous tail;
       cross-checked against the independent zq_L.py at many points of D = -20, -84 up to t = 200, and against the direct lattice sum);
  (ii) an explicit a priori bound M >= |Z_Q| on a disk |w - c| <= R inside Re w > 1/2, from the lattice-point identity (PROVED)
       Z_Q(s) = sum_{n<=X} r(n) n^{-s} + A X^{1-s}/(s-1) - P(X) X^{-s} + s int_X^inf P(u) u^{-s-1} du,   Re s > 1/2,
       with N(u) = #{v != 0: Q(v) <= u} = A u + P(u), A = 2 pi/sqrt|D|, |P(u)| <= 2 A delta sqrt(u) + A delta^2 + 1
       (unit-cell argument, delta = max_{corners} sqrt Q(+-1/2, +-1/2)). So
       |Z_Q(w)| <= sum_{n<=X} r(n) n^{-s_lo} + A X^{1-s_lo}/|w-1|_min + |P(X)| X^{-s_lo} + |w|_max [c1 X^{1/2-s_lo}/(s_lo-1/2) + c0 X^{-s_lo}/s_lo],
       s_lo = Re c - R, c1 = 2 A delta, c0 = A delta^2 + 1. The finite sum is a float64 sum of positive terms, inflated by 1e-8 relative
       (rounding is <= 3e-10 relative for X = 2e6 terms).
Cauchy: Z(w) = Z(c) + a1 (w-c) + T(w), |T| <= M (rho/R)^2/(1-rho/R) on |w-c| = rho. a~ = (Z(c+h)-Z(c-h))/2h has |a1 - a~| <= M (h/R)^2/(R(1-(h/R)^2)).
Rouche against g(w) = Z(c) + a~ (w-c): if |a~| rho - |Z(c)| > |a1-a~| rho + M (rho/R)^2/(1-rho/R), then Z_Q has exactly one zero
(with multiplicity) in |w - c| < rho.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, math, time, numpy as np
sys.path.insert(0,'.')
import eps
from flint import acb, arb, ctx
ctx.prec=96
X=2_000_000
FORMS=[(1,1,6),(2,1,3),(1,0,14),(2,0,7),(3,2,5)]
def setup(Q):
    a,b,c=Q; D=b*b-4*a*c; lam=((a+c)-math.sqrt((a-c)**2+b*b))/2; M=int(math.sqrt(X/lam))+3
    r=np.zeros(X+1,dtype=np.int64)
    for x in range(-M,M+1):
        y=np.arange(-M,M+1); v=a*x*x+b*x*y+c*y*y; v=v[(v>0)&(v<=X)]; np.add.at(r,v,1)
    ns=np.nonzero(r)[0]
    A=2*math.pi/math.sqrt(-D); dl=math.sqrt(max(a+b+c,a-b+c)/4)
    PX=abs(int(r.sum())-A*X)
    return dict(ns=ns.astype(float),rn=r[ns].astype(float),A=A,dl=dl,PX=PX)
def Mbound(S,c,R):
    slo=c.real-R; assert slo>0.5
    part=float(np.sum(S['rn']*S['ns']**(-slo)))*(1+1e-8)
    A,dl=S['A'],S['dl']; c1=2*A*dl; c0=A*dl*dl+1
    d1=abs(c-1)-R; wmax=abs(c)+R
    return part + A*X**(1-slo)/d1 + (S['PX']*(1+1e-9)+1e-6)*X**(-slo) + wmax*(c1*X**(0.5-slo)/(slo-0.5) + c0*X**(-slo)/slo)
def certify(Q,S,x,y,rho=1e-8,h=1e-5):
    c=complex(x,y); R=min(0.1,(x-0.5)/2)
    M=Mbound(S,c,R)
    Zc=eps.Z_Q(Q,acb(x,y)); Zp=eps.Z_Q(Q,acb(x+h,y)); Zm=eps.Z_Q(Q,acb(x-h,y))
    at=(Zp-Zm)/(2*h)
    ea=M*(h/R)**2/(R*(1-(h/R)**2))
    T=M*(rho/R)**2/(1-rho/R)
    lhs=abs(at)*rho-abs(Zc); rhs=arb(ea*rho+T)*(1+arb(1e-12))
    ok=bool(lhs>rhs)
    return dict(Q=list(Q),center=[x,y],rho=rho,R=R,M=M,absZc_upper=float(abs(Zc).upper()),abs_a=float(abs(at).mid()),
                margin=float((lhs-rhs).lower()),certified_one_zero=ok)
if __name__=='__main__':
    zs={}
    for b in json.load(open('scan_T100.json')):
        if tuple(b['Q']) in FORMS: zs.setdefault(tuple(b['Q']),[]).extend(b['zeros'])
    out=[]
    for Q in FORMS:
        t=time.time(); S=setup(Q); print(Q,'setup',round(time.time()-t,1),'s, zeros:',len(zs.get(Q,[])),flush=True)
        for z in zs.get(Q,[]):
            r=certify(Q,S,z[0],z[1]); out.append(r); print(r,flush=True)
    json.dump(out,open('certified_zeros_rouche_D23_D56.json','w'),indent=1)
    print('certified',sum(r['certified_one_zero'] for r in out),'of',len(out))
