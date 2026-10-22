# check sum_{m!=0} Delta_m^{-s} (pure momentum, w=0) = (rho2 sqrt|D|)^s Z_Q(s) at s=3, rho2 in {1, 1.7}
import sys; sys.path.insert(0,'.')
import mpmath as mp
from narain import setup
from theta import Lambda
from flint import acb, arb, ctx
ctx.prec=80
for Q in [(1,0,5),(2,2,3),(1,1,6),(3,2,5),(5,4,5)]:
    for rho2 in [1,1.7]:
        G,Gi,D=setup(Q,rho2)
        M=120; s=3; tot=mp.mpf(0)
        for m1 in range(-M,M+1):
            for m2 in range(-M,M+1):
                if m1==0 and m2==0: continue
                v=mp.matrix([m1,m2]); tot+=((v.T*Gi*v)[0]/2)**(-s)
        L=Lambda(Q,acb(s)); g=(arb(-D).sqrt()/(2*arb.pi()))**3*acb(3).gamma()
        ZQ=L/g
        print(Q,rho2, mp.nstr(tot,12), (arb(rho2)*arb(-D).sqrt())**3*ZQ.real)
