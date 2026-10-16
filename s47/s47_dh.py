"""S47 H47a: ACTUAL Davenport-Heilbronn Weil form on the S14 exact Arb assembly (Legendre*(1-u^2)^4, g in [-L,L]).
Difference trick (no new arch assembly needed):
  arch_DH - arch_zeta = (1/2pi) int Fhat(r)[Re psi(3/4+ir/2) - Re psi(1/4+ir/2) + log 5] dr
                      = log5*F(0) + int_0^{2L} F(x) sech(x/2) dx      [psi(3/4+z)-psi(1/4+z) <-> (1/2)sech(x/2), PROVED below]
  W_DH = (Pole+Arch)_zeta - Pole_zeta + log5*F(0) + int_0^{2L} F sech(x/2) - sum_{2<=n<e^{2L}} 2 b(n)/sqrt(n) F(log n)
  b = (a log) * a^{-1}, a(n) = Re chi(n) + kappa Im chi(n), chi mod 5, chi(2)=i, kappa=(sqrt(10-2sqrt5)-2)/(sqrt5-1)
  (kappa = sqrt(phi+2)-phi: phi from the INPUT CHOICE q=5).
Identity: psi(3/4+ir/2)-psi(1/4+ir/2) = int_0^inf (e^{-t/4}-e^{-3t/4}) e^{-irt/2}/(1-e^{-t}) dt = int_0^inf e^{-irx}/cosh(x/2) dx,
  real part = pi sech(pi r) = FT of (1/2)sech(x/2).  sech moments by acb.integral (rigorous)."""
import sys, json, time; sys.path.insert(0,'.'); A_=sys.argv[:]; sys.argv=sys.argv[:1]
from s46_planted import *
sys.argv=A_
def chi5(n):
    return {0:(0,0),1:(1,0),2:(0,1),3:(0,-1),4:(-1,0)}[n%5]
def bcoef(N):
    s5=arb(5).sqrt(); kap=((10-2*s5).sqrt()-2)/(s5-1)
    a={n:arb(chi5(n)[0])+kap*chi5(n)[1] for n in range(1,N+1)}
    b={1:arb(0)}
    for n in range(2,N+1):
        b[n]=a[n]*arb(n).log()-sum((b[d]*a[n//d] for d in range(2,n) if n%d==0),arb(0))
    return b,kap
def W_DH(M,h,prec=768):
    A,T,ns=build_all(M,h,prec); ctx.prec=prec; ha=arb(h); hc=acb(ha)
    Gs=Gpolys(M); D=max(g.degree() for g in Gs.values())+1
    ms=[acb.integral(lambda z,an,m=m: z**m/((hc*z/2).cosh()),0,2).real for m in range(D+1)]
    m3=[acb.integral(lambda z,an,m=m: z**m*(hc*z/2).cosh(),0,2).real for m in range(D+1)]  # zeta pole moments (as S14)
    N=int(mp.floor(mp.exp(2*mp.mpf(h)))); nn=[n for n in range(2,N+1) if arb(n).log()<2*ha]
    b,kap=bcoef(max(nn)); l5=arb(5).log()
    W=arb_mat(M,M)
    for (j,k),g in Gs.items():
        co=[arb(c.p)/c.q for c in g.coeffs()]
        sech=ha*ha*sum(co[m]*ms[m] for m in range(len(co)))
        pr=arb(0)
        for n in nn:
            v=arb(n).log()/ha; pr+=2*b[n]/arb(n).sqrt()*ha*sum(co[m]*v**m for m in range(len(co)))
        pole=4*ha*ha*sum(co[m]*m3[m] for m in range(len(co)))
        W[j,k]=A[j,k]-pole+l5*ha*co[0]+sech-pr
    Wz=arb_mat(M,M)
    for i in range(M):
        for j in range(M):
            Wz[i,j]=A[i,j]-sum((2*arb(vonmangoldt(n)).log()/arb(n).sqrt())*T[n][i,j] for n in ns)
    return W,Wz,nn,b,kap
if __name__=="__main__":
    M=int(sys.argv[1]); hs=sys.argv[2].split(','); mp.mp.dps=int(sys.argv[3]) if len(sys.argv)>3 else 100
    prec=int(sys.argv[4]) if len(sys.argv)>4 else 768
    for h in hs:
        t0=time.time(); W,Wz,nn,b,kap=W_DH(M,h,prec)
        row={"M":M,"L":h,"kind":"DH_actual","n_used":nn,"kappa":kap.str(15),
             "maxrad":max(float(W[i,j].rad()) for i in range(M) for j in range(M))}
        for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
            lz,_=secmin(Wz,I); row["zeta_"+nm]=mp.nstr(lz,6)
            l,v=secmin(W,I); row["DH_"+nm]=mp.nstr(l,8)
            if l<0:
                val=certify(W,I,v); row[nm+"_cert"]=val.str(8); row[nm+"_certified_neg"]=bool(val<0)
        row["sec"]=round(time.time()-t0,1); print(json.dumps(row),flush=True)
