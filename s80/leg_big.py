import time, math, json, sys
from flint import arb, ctx
sys.argv=[sys.argv[0]]
exec(open('s80/leg_fixedprec.py').read().split('N=int')[0])
N=20000; res={}
for yv in ['0.5','0.99','0.9999','0.9999999']:
    y=arb(yv); fy=float(yv); boost=int(N*math.log2(fy+math.sqrt(1+fy*fy)))+64+320
    t=time.time(); ref=fwd(y,N,boost+200); tr=time.time()-t
    t=time.time(); P=fwd(y,N,320); tp=time.time()-t
    t=time.time(); fwd(y,N,boost); tb=time.time()-t
    err=max(abs(float((P[d]-ref[d]).mid())) for d in range(N)); u=2.0**-320
    s=math.sqrt(1-fy*fy); res[yv]=dict(t320=tp,tboost=tb,boostbits=boost,err_over_u=err/u,bound_guess=N/s)
    print(yv,'t320=%.3fs tboost(%d bits)=%.3fs speedup %.1fx  maxerr/u=%.3g  N/sin=%.3g'%(tp,boost,tb,tb/tp,err/u,N/s),flush=True)
json.dump(res,open('s80/leg_N20000.json','w'),indent=1)
