"""S16: n=3 window, ½log3<h<log2. W = A - c2 T2 - c3 T3 on C_h (S14 exact Arb basis, Legendre*(1-u^2)^4).
Even/odd sectors decouple. Admissible region R_M(h) = {(c2,c3): A-c2T2-c3T3 PSD} (Ritz, convex).
Boundary computed radially from zeta point: B=W_zeta (PD on Ritz space), edge along d = 1/lambda_max(B^{-1/2}(d.T)B^{-1/2})."""
import sys, json; sys.path.insert(0,'../s14')
from s14_exact import *; import mpmath as mp
mp.mp.dps=80
def build3(M,h):
    A,T2=build(M,h,512); ctx.prec=512; hh=arb(h); v3=arb(3).log()/hh
    T3=arb_mat(M,M)
    for j in range(M):
        for k in range(j,M):
            if (j+k)%2==0:
                g=(G(j,k)+G(k,j))/2; co=[arb(c.p)/c.q for c in g.coeffs()]
                T3[j,k]=T3[k,j]=hh*sum(co[m]*v3**m for m in range(len(co)))
    return A,T2,T3
def tomp(Ma,I): return mp.matrix([[mp.mpf(Ma[i,j].mid().str(100,radius=False)) for j in I] for i in I])
M=int(sys.argv[1]); h=sys.argv[2]
A,T2,T3=build3(M,h)
c2z=mp.sqrt(2)*mp.log(2); c3z=2*mp.log(3)/mp.sqrt(3)
th=mp.mpf(1)/4; c2F=c2z*mp.cosh(th*mp.log(2)); c3F=c3z*mp.cosh(th*mp.log(3))
res={"M":M,"h":h,"zeta":[str(c2z),str(c3z)],"F14":[mp.nstr(c2F,12),mp.nstr(c3F,12)]}
for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
    a,t2,t3=tomp(A,I),tomp(T2,I),tomp(T3,I); B=a-c2z*t2-c3z*t3
    ev,V=mp.eigsy(B); 
    if min(ev)<=0: res[nm]="zeta not PD"; continue
    S=V*mp.diag([1/mp.sqrt(e) for e in ev])*V.T
    X2=S*t2*S; X3=S*t3*S
    def edge(d2,d3):
        e=mp.eigsy(d2*X2+d3*X3,eigvals_only=True); lm=max(e)
        return (1/lm) if lm>0 else mp.inf
    # axis edges & F point test
    r={"c2_up":mp.nstr(edge(1,0),8),"c2_dn":mp.nstr(-edge(-1,0),8),"c3_up":mp.nstr(edge(0,1),8),"c3_dn":mp.nstr(-edge(0,-1),8),
       "minEig_B":mp.nstr(min(ev),5)}
    pts=[]
    for k in range(72):
        ang=2*mp.pi*k/72; d2,d3=mp.cos(ang),mp.sin(ang); t=edge(d2,d3)
        pts.append([float(ang),float(t) if t!=mp.inf else None])
    r["radial"]=pts
    dF=(c2F-c2z,c3F-c3z); tF=edge(*dF); r["F14_inside"]= bool(tF>1); r["F14_edge_frac"]=mp.nstr(tF,8)
    # rank of T2,T3 in sector
    r["eigT2"]=[mp.nstr(x,4) for x in mp.eigsy(t2,eigvals_only=True) if abs(x)>1e-30]
    r["eigT3"]=[mp.nstr(x,4) for x in mp.eigsy(t3,eigvals_only=True) if abs(x)>1e-30]
    res[nm]=r
print(json.dumps(res))
