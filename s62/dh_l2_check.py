# cross-check of the cosritz FORM=DH assembly against the independent S47 Legendre*(1-x^2)^4 DH assembly, L2-normalised Ritz (NUMERICAL)
import sys; sys.path.insert(0,'.'); ARGS=sys.argv[:]; sys.argv=sys.argv[:1]
from s47_dh import *
M=int(ARGS[1]); h=ARGS[2]; mp.mp.dps=int(ARGS[3]); prec=int(ARGS[4])
W,Wz,nn,b,kap=W_DH(M,h,prec); ctx.prec=prec; ha=arb(h); Gs=Gpolys(M)
Gr=arb_mat(M,M)
for i in range(M):
    for j in range(M):
        if (i+j)%2==0: c=Gs[i,j].coeffs(); Gr[i,j]=ha*(arb(c[0].p)/c[0].q)
out={}
for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
    for tag,X in (("DH",W),("zeta",Wz)):
        Wm=tomp(X,I); G=tomp(Gr,I); Lc=mp.cholesky(G); Li=mp.inverse(Lc); C=Li*Wm*Li.T; C=(C+C.T)/2
        out[tag+"_"+nm]=mp.nstr(min(mp.eigsy(C,eigvals_only=True)),8)
print(M,h,out)
