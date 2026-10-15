"""S16 floor test: which piece supplies the odd-sector floor / even cap. For each sector compute edges of
c2 (c3 fixed at zeta) for W = s_pole*Pole + Arch - c2T2 - c3T3, scanning s_pole. Also M-convergence of edges."""
import sys,json; sys.path.insert(0,'.')
from s16_exactP import *; import mpmath as mp
mp.mp.dps=80
M=int(sys.argv[1]); h=sys.argv[2]
A,T2,P=build(M,h,512); ctx.prec=512; hh=arb(h); v3=arb(3).log()/hh; T3=arb_mat(M,M)
for j in range(M):
    for k in range(j,M):
        if (j+k)%2==0:
            g=(G(j,k)+G(k,j))/2; co=[arb(c.p)/c.q for c in g.coeffs()]
            T3[j,k]=T3[k,j]=hh*sum(co[m]*v3**m for m in range(len(co)))
tm=lambda X,I: mp.matrix([[mp.mpf(X[i,j].mid().str(100,radius=False)) for j in I] for i in I])
c2z=mp.sqrt(2)*mp.log(2); c3z=2*mp.log(3)/mp.sqrt(3)
def edges(B,t):
    ev,V=mp.eigsy(B)
    if min(ev)<=0: return None
    S=V*mp.diag([1/mp.sqrt(e) for e in ev])*V.T; e=mp.eigsy(S*t*S,eigvals_only=True)
    return (mp.nstr(1/max(e),6) if max(e)>0 else 'inf', mp.nstr(1/min(e),6) if min(e)<0 else '-inf')
out={"M":M,"h":h}
for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
    a,p,t2,t3=tm(A,I),tm(P,I),tm(T2,I),tm(T3,I); r={}
    for sp in ('1','1.0001','0.9999','1.001','0.999'):
        sp=mp.mpf(sp); B=a+(sp-1)*p-c2z*t2-c3z*t3
        r[mp.nstr(sp,5)]=edges(B,t2)  # (up,down) offsets of c2 from zeta value
    r["minEig_W"]=mp.nstr(min(mp.eigsy(a-c2z*t2-c3z*t3,eigvals_only=True)),5)
    r["minEig_arch+pole"]=mp.nstr(min(mp.eigsy(a,eigvals_only=True)),5)
    out[nm]=r
print(json.dumps(out))
# kappa window: W(kappa)=A+(kappa-1)P - c2z T2 - c3z T3 ; edges in kappa per sector
kw={}
for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
    a,p,t2,t3=tm(A,I),tm(P,I),tm(T2,I),tm(T3,I)
    e=edges(a-c2z*t2-c3z*t3,-p)   # B - t*(-P) = B + tP
    kw[nm]=e; kw[nm+"_Psign"]=[mp.nstr(min(mp.eigsy(p,eigvals_only=True)),4),mp.nstr(max(mp.eigsy(p,eigvals_only=True)),4)]
print(json.dumps({"kappa_minus_1_edges(up,down)":kw}))
