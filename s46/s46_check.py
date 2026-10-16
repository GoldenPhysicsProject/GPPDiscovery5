"""Convention check (NUMERICAL): sum over first N zeta zeros of 2 int F cos(gamma x) dx vs W_zeta on the S14 assembly."""
import sys; sys.path.insert(0,'.')
from s46_planted import *
M=int(sys.argv[1]); h=sys.argv[2]; N=int(sys.argv[3]); prec=512
mp.mp.dps=40
W,ns=Wzeta(M,h,prec); Gs=Gpolys(M)
Z=zeros(N)
S=arb_mat(M,M)
ctx.prec=512
D=max(g.degree() for g in Gs.values())+1
mu=[arb(0)]*(D+1); ha=arb(h)
for z in Z:
    I=moms(acb(0,1)*ha*arb(z),D)
    mu=[mu[m]+2*I[m].real for m in range(D+1)]
for (j,k),g in Gs.items():
    co=[arb(c.p)/c.q for c in g.coeffs()]
    S[j,k]=2*ha*ha*sum(co[m]*mu[m] for m in range(len(co)))
mx=0
for i in range(M):
    for j in range(M):
        if (i+j)%2==0:
            w=float(W[i,j].mid()); s=float(S[i,j].mid()); mx=max(mx,abs(w-s)/max(abs(w),1e-300))
            if i<3 and j<3: print(i,j,w,s)
print("max rel diff",mx)
