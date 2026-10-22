"""NUMERICAL check of the c=2 Narain scalar-primary spectral zeta identity
   S(s) := sum'_{scalar U(1)^2 primaries} Delta^{-s} = |D|^{s/2} * E_s(rho) * Z_Q(s)     (B=0, rho = i*rho2, sigma = tau_Q)
It follows from Benjamin-Collier-Fitzpatrick-Maloney-Perlmutter arXiv:2107.10744 eq. (3.25),
[sum' delta_{n.w,0} (M^2)^{-s} = 2 zeta(2s) E_s(rho) E_s(sigma), M^2 = 2 Delta],
together with 2 zeta(2s) E_s(tau_Q) = (sqrt|D|/2)^s Z_Q(s), where E_s is the coprime Poincare sum.
Here 2 zeta(2s) E_s(i*y) = y^s * sum'_{c,d} (d^2 + y^2 c^2)^{-s}  (computed by direct summation at s=3).
Scalars (J = m.w = 0) are enumerated as: pure momentum; pure winding; and m = d*u, w = c*(-u2,u1)
with u primitive, c >= 1, d != 0 (the BCFMP parametrisation, eq. 3.23).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys; sys.path.insert(0,'.')
from math import gcd
import numpy as np, mpmath as mp
from narain import setup
from zq_L import ZQ
from theta import Lambda
from flint import acb, arb, ctx
ctx.prec=90
def ZQ_any(Q,s):
    a,b,c=Q; D=b*b-4*a*c
    g=(arb(-D).sqrt()/(2*arb.pi()))**s*acb(s).gamma()
    return (Lambda(Q,acb(s),bits=80)/g).real
s=3.0; M=150
out=[]
for Q in [(1,0,5),(2,2,3),(1,1,6),(3,2,5)]:
    for rho2 in [1.0,1.3]:
        G,Gi,D=setup(Q,rho2)
        Gf=np.array([[float(G[i,j]) for j in range(2)] for i in range(2)]); Gif=np.linalg.inv(Gf)
        r=np.arange(-M,M+1); a1,a2=np.meshgrid(r,r,indexing='ij'); a1=a1.ravel().astype(float); a2=a2.ravel().astype(float)
        nz=(a1!=0)|(a2!=0)
        mom=0.5*(Gif[0,0]*a1**2+2*Gif[0,1]*a1*a2+Gif[1,1]*a2**2)[nz]
        win=0.5*(Gf[0,0]*a1**2+2*Gf[0,1]*a1*a2+Gf[1,1]*a2**2)[nz]
        tot=np.sum(mom**-s)+np.sum(win**-s)
        prim=np.array([np.gcd(int(x),int(y))==1 for x,y in zip(a1[nz],a2[nz])])
        mu=mom[prim]; u1=a1[nz][prim]; u2=a2[nz][prim]
        wu=0.5*(Gf[0,0]*u2**2-2*Gf[0,1]*u1*u2+Gf[1,1]*u1**2)
        K=40
        cc=np.arange(1,K+1,dtype=float); dd=np.concatenate([np.arange(-K,0),np.arange(1,K+1)]).astype(float)
        CC,DD=np.meshgrid(cc,dd,indexing='ij'); CC=CC.ravel(); DD=DD.ravel()
        for i in range(len(mu)):
            tot+=np.sum((DD*DD*mu[i]+CC*CC*wu[i])**-s)
        L2=400; rr=np.arange(-L2,L2+1,dtype=float); X,Y=np.meshgrid(rr,rr); V=(Y**2+rho2**2*X**2).ravel(); V=V[V>0]
        Es_rho = arb(rho2**3*float(np.sum(V**-3.0)))/(2*arb(6).zeta())   # truncation ~1e-9 relative
        rhs = arb(-D)**arb(1.5)*Es_rho*ZQ_any(Q,3)
        line=f"{Q} rho2={rho2}  lattice S(3)={tot:.10g}   |D|^(3/2) E_3(rho) Z_Q(3)={rhs.str(10)}"
        print(line,flush=True); out.append(line)
open('scalar_zeta_identity.log','w').write("\n".join(out)+"\n")
