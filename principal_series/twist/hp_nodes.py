"""F_theta failing even vector node count vs N (mpmath dps 80, S87 builder). HIGH-PRECISION NUMERICAL. RH OPEN.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json; sys.argv=['x','24']; sys.path.insert(0,'/workspace/newbot/twist')
import hp_onset_lib as H, mpmath as mp, numpy as np
from numpy.polynomial import legendre as leg
res=[]
xs=np.linspace(0,1,40001)
for L in ['0.5','0.6','0.7','0.8']:
    for N in [24,32,40]:
        Qe,Qo=H.forms('FT', mp.mpf(L), N, mp.mpf('0.1'))
        E,V=mp.eigsy(mp.matrix(Qe)); k=min(range(N),key=lambda i:E[i])
        v=np.array([float(V[i,k]) for i in range(N)])
        f=sum(v[i]*np.sqrt(2*i+1)*leg.legval(2*xs-1,[0]*i+[1]) for i in range(N))
        g=f[np.abs(f)>1e-10*np.abs(f).max()]; nodes=int(((g[:-1]*g[1:])<0).sum())
        r={"theta":0.1,"L":L,"N":N,"lam_min_even":mp.nstr(E[k],8),"nodes_on_[0,L]":nodes}
        print(r,flush=True); res.append(r)
json.dump(res,open('hp_nodes_F.json','w'),indent=1)
