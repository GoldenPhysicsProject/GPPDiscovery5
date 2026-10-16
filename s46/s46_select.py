"""Selected-config scan: argv M hs_csv cfgname gamma_or_mean ia ib dels_csv [dps]. Same pipeline as s46_scan (Arb W + Arb Delta, mp eig, Arb witness)."""
import sys; sys.path.insert(0,'.'); A_=sys.argv[:]; sys.argv=sys.argv[:1]
from s46_planted import *
sys.argv=A_
M=int(sys.argv[1]); hs=sys.argv[2].split(','); name=sys.argv[3]; gp=sys.argv[4]; ia,ib=int(sys.argv[5]),int(sys.argv[6]); dels=sys.argv[7].split(',')
mp.mp.dps=int(sys.argv[8]) if len(sys.argv)>8 else 100
Z=zeros(60); Gs=Gpolys(M)
for h in hs:
    t0=time.time(); W,ns=Wzeta(M,h,768)
    ga,gb=zball(Z[ia-1]),zball(Z[ib-1]); gam=(ga+gb)/2 if gp=='mean' else arb(gp)
    row={"M":M,"h":h,"cfg":name}
    for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
        l,_=secmin(W,I); row["zeta_"+nm]=mp.nstr(l,6)
    for d in dels:
        X=W+delta_mat(M,h,Gs,planted_terms(gam,arb(d),ga,gb)); r={}
        for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
            l,v=secmin(X,I); r[nm]=mp.nstr(l,6)
            if l<0: val=certify(X,I,v); r[nm+"_cert"]=val.str(8); r[nm+"_certified_neg"]=bool(val<0)
        row[d]=r
    row["sec"]=round(time.time()-t0,1); print(json.dumps(row),flush=True)
