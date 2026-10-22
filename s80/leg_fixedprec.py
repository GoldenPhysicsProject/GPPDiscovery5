"""Grok S80 H80a: forward Legendre recurrence at FIXED precision p (no ball boost) vs boosted reference.
NUMERICAL test of error growth in degree d, to replace the +N*log2(|y|+sqrt(1+y^2)) bit boost by an a-priori bound."""
import sys, time, math, json
from flint import arb, ctx
def fwd(y, N, prec):
    old=ctx.prec; ctx.prec=prec
    try:
        y=arb(y.mid()); P=[arb(1), y]
        for d in range(1,N-1):
            P.append((((2*d+1)*y*P[d]-d*P[d-1])/(d+1)).mid())  # point arithmetic, discard radius
            P[-1]=arb(P[-1])
    finally: ctx.prec=old
    return P
N=int(sys.argv[1]); out={}
for yv in ['0.3','0.9','0.999','0.99999']:
    ctx.prec=600; y=arb(yv)
    t=time.time(); ref=fwd(y,N,N*2+600); tr=time.time()-t
    for p in (192,320):
        t=time.time(); P=fwd(y,N,p); tp=time.time()-t
        err=[abs(float((P[d]-ref[d]).mid())) for d in range(N)]
        # growth: max err in [d/2,d] relative to 2^-p
        u=2.0**-p
        g=[max(err[max(1,d//2):d+1])/u for d in (100,1000,N-1) if d<N]
        out[f'{yv}_p{p}']=dict(t=tp,tref=tr,growth_over_u=g,maxerr=max(err))
        print(yv,p,'t=%.2fs ref=%.2fs'%(tp,tr),'err/u at d=100,1000,N:',['%.3g'%x for x in g],flush=True)
json.dump(out,open(f's80/leg_N{N}.json','w'),indent=1)
