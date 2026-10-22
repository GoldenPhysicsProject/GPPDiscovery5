import sys; sys.path.insert(0,'.')
from theta import Lambda
from eps import Lambda_Q as LamEis, tau_of
from flint import acb, arb, ctx
import sympy
ctx.prec=80
def L(d,s):
    q=abs(d); s=acb(s); r=acb(0)
    for a in range(1,q+1):
        k=int(sympy.kronecker_symbol(d,a))
        if k: r+=k*s.zeta(acb(a)/q)
    return r*acb(q)**(-s)
def gam(D,s): return (arb(-D).sqrt()/(2*arb.pi()))**s*s.gamma()
for Q in [(1,0,5),(2,2,3),(1,1,6),(2,1,3)]:
    print(Q,'s=3', Lambda(Q,acb(3)), LamEis(Q,acb(3)))
for s in [acb(0.7,30.3), acb(1.3,85), acb(0.62,140.1)]:
    z1=Lambda((1,0,5),s); z2=Lambda((2,2,3),s); g=gam(-20,s)
    print(s,'sum ', z1+z2, 2*g*s.zeta()*L(-20,s))
    print(s,'diff', z1-z2, 2*g*L(-4,s)*L(5,s))
    print('FE', Lambda((1,0,5),1-s))
