"""S14: exact-rational / Arb construction of the n=2-window Weil form A - c T2 on C_h (Legendre*(1-u^2)^4 basis, as S13).
F_jk(x)=h G_jk(x/h), G exact fmpq_poly on [0,2]. Arch term in real space (Gauss digamma integral):
 arch(F) = (-gamma-log pi)F(0) - F(0)log(1-e^{-4h}) + int_0^{2h} [(F0-F)/x] e^{-x} x/sinh x dx - 2 int_0^{2h} F e^{-x/4} sinh(3x/4)/sinh x dx
 pole(F) = int F(x) 2cosh(x/2) dx;  T = F(log 2);  W_zeta = pole+arch - sqrt2 log2 * T.
All moments int_0^2 v^m phi(v) dv by acb.integral (rigorous)."""
from flint import fmpq_poly, fmpq, arb, acb, ctx, arb_mat
import sys, functools
def legendre(n):
    P0,P1=fmpq_poly([1]),fmpq_poly([0,1])
    if n==0: return P0
    for k in range(1,n): P0,P1=P1,(fmpq(2*k+1,k+1)*fmpq_poly([0,1])*P1-fmpq(k,k+1)*P0)
    return P1
import os; EDGEK=int(os.environ.get("EDGEK","4")); W4=fmpq_poly([1,0,-1])**EDGEK
def p(j): return legendre(j)*W4
def G(j,k):
    pj,pk=p(j),p(k); res=fmpq_poly([0]); dk=pk; v=fmpq_poly([0,1]); fact=1
    for i in range(pk.degree()+1):
        Q=(pj*dk).integral()
        res+= (fmpq_poly([0,-1])**i)*fmpq(1,fact)*(Q(fmpq(1)) - Q(v-1))
        dk=dk.derivative(); fact*=(i+1)
    return res
def build(M,h,prec=512):
    ctx.prec=prec; h=arb(h) if not isinstance(h,arb) else h
    Gs={}
    for j in range(M):
        for k in range(j,M):
            if (j+k)%2==0: g=(G(j,k)+G(k,j))/2; Gs[j,k]=Gs[k,j]=g
    D=max(g.degree() for g in Gs.values())+1
    def mom(phi):
        out=[]
        for m in range(D+1):
            f=lambda z,an,m=m: z**m*phi(z)
            out.append(acb.integral(f,0,2).real)
        return out
    hc=acb(h)
    m1=mom(lambda z: (-hc*z).exp()/ (acb(0,1)*hc*z).sinc() )   # e^{-hv} hv/sinh(hv)
    m2=mom(lambda z: (-hc*z/4).exp()*(acb(0,1)*hc*z*3/4).sinc()*fmpq(3,4)/(acb(0,1)*hc*z).sinc())  # e^{-hv/4}sinh(3hv/4)/sinh(hv)
    m3=mom(lambda z: (hc*z/2).cosh())
    L2=arb(2).log(); vl=L2/h
    A=arb_mat(M,M); T=arb_mat(M,M)
    c0=-arb.const_euler()-arb.pi().log()-(1-(-4*h).exp()).log()
    for (j,k),g in Gs.items():
        co=[arb(c.p)/c.q for c in g.coeffs()]
        G0=co[0]; q=[-c for c in co[1:]]
        pole=4*h*h*sum(co[m]*m3[m] for m in range(len(co)))
        arch=h*c0*G0 + h*sum(q[m]*m1[m] for m in range(len(q))) - 2*h*h*sum(co[m]*m2[m] for m in range(len(co)))
        A[j,k]=pole+arch
        T[j,k]=h*sum(co[m]*vl**m for m in range(len(co))) if vl<2 else arb(0)
    return A,T
