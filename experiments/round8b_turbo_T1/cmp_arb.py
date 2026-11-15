import json, glob, numpy as np
from flint import arb
E = {}
for f in sorted(glob.glob('/workspace/r6/GPPDiscovery5/experiments/round7_turbo/shard_*.json')): E.update(json.load(open(f)))
l = arb(1)/4; pi = arb.pi(); c2 = (2/l).sqrt()
om = lambda j: j*pi/l; w = lambda j: c2*(-1)**(j+1)/om(j)
def ent(a,b):
    if a=='d' and b=='d': return arb(E['d,d'])
    if b=='d': return arb(E[f"{a-1},d"])
    if a=='d': return arb(E[f"{b-1},d"])
    a,b=min(a,b),max(a,b); return arb(E[f"{a-1},{b-1}"])
g = ent('d','d')
for j,k in [(13,13),(13,14),(13,20),(30,30),(30,31)]:
    v = ent(j,k)-w(k)*ent(j,'d')-w(j)*ent(k,'d')+w(j)*w(k)*g
    print(j,k,v.str(12))
