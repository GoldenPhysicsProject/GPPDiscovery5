"""Queen S1: sign-change count of the 2nd even and 1st odd eigenvectors of the non-pole operator A-P (zeta, cache s91). NUMERICAL dps 90."""
import pickle, mpmath as mp, numpy as np
from numpy.polynomial import legendre as Lg
mp.mp.dps=90
cache=pickle.load(open("../../s91/forms_cache.pkl","rb"))
def M(c,k): return mp.matrix([[mp.mpf(x) for x in r] for r in c[k]])
def nodes(vec,L,N):
    s=np.linspace(0,L,4001)[1:-1]; t=2*s/L-1
    f=sum(float(vec[j])*np.sqrt((2*j+1)/L)*Lg.legval(t,np.eye(N)[j]) for j in range(N))
    sg=np.sign(f[np.abs(f)>1e-6*np.abs(f).max()]); return int(np.sum(sg[1:]!=sg[:-1]))
for key in sorted(cache):
    F,th,L,N=key
    if F!="ZETA": continue
    c=cache[key]; row=[L,N]
    for par in "eo":
        E,V=mp.eigsy(M(c,"A"+par)-M(c,"P"+par)); idx=sorted(range(N),key=lambda i:E[i])
        row+=[par]+[nodes([V[j,idx[r]] for j in range(N)],L,N) for r in range(3)]
    print(row,flush=True)
