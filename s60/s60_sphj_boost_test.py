# Grok S60: diagnose S58 even L=7/5 T=1310 failure: precision of sph_j_all backward recurrence at x = T*L (or nearby), D ~ 2994
import os, sys, time
sys.path.insert(0, '/workspace/GPPDiscovery5/s53'); sys.path.insert(0, '/workspace/GPPDiscovery5/s27'); sys.path.insert(0, '/workspace/GPPDiscovery5/s36')
from flint import arb, ctx
ctx.prec = 384
import s53_onestroke_par as S
for xb in [1834, 1310*1.4*0.9]:
    for B in [1.75, 2.0, 2.5]:
        S.BOOST = B; t0 = time.time()
        v = S.sph_j_all(arb(xb), 2994)
        rads = [float(t.rad()) for t in v]
        print(f'x={xb:.1f} BOOST={B}: max rad={max(rads):.3e} rad j0={rads[0]:.3e} j0={v[0].str(5)} [{time.time()-t0:.1f}s]', flush=True)
