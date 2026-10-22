"""Rigorous ball evaluation of Z_Q(s) for forms whose class group is 2-torsion (only genus characters),
via the classical decomposition (Dirichlet 1840 / Kronecker; see e.g. Y. Lee, arXiv:1204.6297, eq. (1.2)):
    Z_Q(s) = (w/h) * sum_{genus chars chi} chi(Q) L(s,chi_{d1}) L(s,chi_{d2}),   D = d1*d2 (fundamental).
L(s,chi_d) = |d|^{-s} sum_a chi_d(a) zeta(s, a/|d|) with Arb Hurwitz zeta (rigorous for ball s).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
from flint import acb, arb, ctx
from sympy import kronecker_symbol
from math import gcd

GENUS = {  # D: list of (d1,d2) fundamental factorisations (unordered), first is (1,D)
 -4:  [(1,-4)],
 -20: [(1,-20),(-4,5)],
 -84: [(1,-84),(-4,21),(-3,28),(-7,12)],
}
W = {-4:4,-20:2,-84:2}
H = {-4:1,-20:2,-84:4}

def Ld(d,s):
    if d==1: return s.zeta()
    q=abs(d); tot=acb(0)
    for a in range(1,q):
        c=int(kronecker_symbol(d,a))
        if c: tot+= c*s.zeta(acb(a)/q)
    return tot*acb(q)**(-s)

def chi_on_form(d1,d2,Q):
    a,b,c=Q
    for (x,y) in [(1,0),(0,1),(1,1),(1,-1),(2,1),(1,2),(3,1),(1,3)]:
        n=a*x*x+b*x*y+c*y*y
        for d in (d1,d2):
            if d!=1 and gcd(n,abs(d))==1:
                return int(kronecker_symbol(d,n))
    raise ValueError

def ZQ(Q,s):
    a,b,c=Q; D=b*b-4*a*c
    tot=acb(0)
    for (d1,d2) in GENUS[D]:
        chi=1 if d1==1 else chi_on_form(d1,d2,Q)
        tot+= chi*Ld(d1,s)*Ld(d2,s)
    return tot*acb(W[D])/H[D]
