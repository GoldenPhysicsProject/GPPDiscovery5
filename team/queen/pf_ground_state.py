import pickle, mpmath as mp, numpy as np, sys
from numpy.polynomial import legendre as Lg
mp.mp.dps=90
cache=pickle.load(open("../../s91/forms_cache.pkl","rb"))
print(sorted(set((k[0],k[1]) for k in cache)), sorted(set(k[2] for k in cache)))
def M(c,k): return mp.matrix([[mp.mpf(x) for x in r] for r in c[k]])
def nodes(vec,L,N):
    s=np.linspace(0,L,4001)[1:-1]; t=2*s/L-1
    f=sum(float(vec[j])*np.sqrt((2*j+1)/L)*Lg.legval(t,np.eye(N)[j]) for j in range(N))
    sg=np.sign(f[np.abs(f)>1e-8*np.abs(f).max()]); return int(np.sum(sg[1:]!=sg[:-1])), float(np.sign(f[len(f)//2]))
out=[]
for key in sorted(cache):
    F,th,L,N=key
    c=cache[key]
    row=[F,th,L,N]
    for par in "eo":
        A=M(c,"A"+par); P=M(c,"P"+par); R=A-P
        E,V=mp.eigsy(R)
        ev=[E[i] for i in range(N)]; idx=sorted(range(N),key=lambda i:ev[i])
        lmax=max(abs(x) for x in ev)
        nneg=sum(1 for x in ev if x<0)
        v0=[V[j,idx[0]] for j in range(N)]
        nd=nodes(v0,L,N)
        pv=c["p"+par]; npv=len(pv)
        # overlap of ground state with pole vector
        ov=None
        if npv:
            p=mp.matrix([mp.mpf(x) for x in pv[0]]); ov=float(sum(p[j]*v0[j] for j in range(N))/mp.sqrt(sum(p[j]**2 for j in range(N))))
        row+= [par, nneg, mp.nstr(ev[idx[0]]/lmax,4), mp.nstr(ev[idx[1]]/lmax,4), nd[0], npv, ov]
    print(row, flush=True)
