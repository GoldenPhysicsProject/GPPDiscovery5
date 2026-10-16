"""S46 H46a scan: lambda_min per sector of Q_fake=W_zeta+Delta on L grid; Arb witness at first detected L per config.
argv: M hs_csv [dps] [prec]"""
import sys; sys.path.insert(0,'.')
from s46_planted import *
M=int(sys.argv[1]); hs=sys.argv[2].split(','); mp.mp.dps=int(sys.argv[3]) if len(sys.argv)>3 else 80
prec=int(sys.argv[4]) if len(sys.argv)>4 else 768
Z=zeros(60)
# configs: (name, gamma_planted, (idx_a, idx_b)) 1-based zero indices
cfg=[("pair1-2",None,(1,2)),("pair10-11",None,(10,11)),("pair50-51",None,(50,51)),("DH85.699",'85.699',(23,24))]
dels=['0.05','0.1','0.2','0.308','0.4']
Gs=Gpolys(M); done=set()
for h in hs:
    t0=time.time(); W,ns=Wzeta(M,h,prec)
    row={"M":M,"h":h,"n_prime_powers":len(ns)}
    for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
        l,_=secmin(W,I); row["zeta_"+nm]=mp.nstr(l,6)
    res={}
    for name,gp,(ia,ib) in cfg:
        ga,gb=zball(Z[ia-1]),zball(Z[ib-1])
        gam=arb(gp) if gp else (ga+gb)/2
        for d in dels:
            X=W+delta_mat(M,h,Gs,planted_terms(gam,arb(d),ga,gb))
            r={}
            for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
                l,v=secmin(X,I); r[nm]=mp.nstr(l,6)
                key=(name,d,nm)
                if l<0 and key not in done:
                    val=certify(X,I,v); r[nm+"_cert"]=val.str(8); r[nm+"_certified_neg"]=bool(val<0)
                    if val<0: done.add(key)
            res[f"{name}|{d}"]=r
    row["configs"]=res; row["sec"]=round(time.time()-t0,1)
    print(json.dumps(row),flush=True)
