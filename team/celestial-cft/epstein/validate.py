import sys; sys.path.insert(0,'.')
from eps import *
from flint import acb, arb, ctx
ctx.prec=80
# 1. against direct lattice sum at s=3
for Q in [(1,0,1),(1,0,5),(2,2,3),(2,1,3),(1,1,6)]:
    print(Q, Z_Q(Q,3), Z_direct(Q,3.0,300))
# 2. against genus/Dirichlet decomposition for D=-20:
#   Z_[1,0,5]+Z_[2,2,3] = 2 zeta L(chi_-20), Z_[1,0,5]-Z_[2,2,3] = 2 L(chi_-4) L(chi_5)
def kron(d,n):
    from sympy import jacobi_symbol
    from sympy.ntheory import legendre_symbol
    from sympy.ntheory.residue_ntheory import jacobi_symbol as js
    import sympy
    return sympy.kronecker_symbol(d,n)
def L(d,s):
    q=abs(d); s=acb(s); r=acb(0)
    for a in range(1,q+1):
        k=kron(d,a)
        if k: r+=int(k)*s.zeta(acb(a)/q)
    return r*acb(q)**(-s)
for s in [acb(0.7,30.3), acb(0.5,14.2), acb(1.3,85)]:
    z1=Z_Q((1,0,5),s); z2=Z_Q((2,2,3),s)
    print('D=-20 sum ', z1+z2, 2*s.zeta()*L(-20,s))
    print('D=-20 diff', z1-z2, 2*L(-4,s)*L(5,s))
    print('FE check', Lambda_Q((1,0,5),s), Lambda_Q((1,0,5),1-s))
