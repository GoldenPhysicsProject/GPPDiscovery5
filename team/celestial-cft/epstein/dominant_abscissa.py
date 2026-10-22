"""Dominant-term abscissa bound. CERTIFIED (Arb ball sums + rigorous tail at the two reported endpoints).
For a positive-definite form Q with minimum m and multiplicity r(m):
  Z_Q(s) = r(m) m^{-s} (1 + E(s)),  |E(s)| <= B(sigma) := sum_{n>m} (r(n)/r(m)) (m/n)^sigma,
B is strictly decreasing, so with sigma_dom := the unique root of B(sigma) = 1 the half-plane Re s > sigma_dom is zero-free,
hence Theta(Q) := sup Re rho <= sigma_dom. We certify sigma_up with B(sigma_up) < 1 (rigorous upper bound incl. tail) and
sigma_lo with B(sigma_lo) > 1 (rigorous lower bound from the finite partial sum), so sigma_dom lies in [sigma_lo, sigma_up].
B(sigma) = 1 is only a sufficient-condition edge: Theta can sit well below sigma_dom.
Tail: put the unit cell v+[-1/2,1/2]^2 around each lattice point; all cells of points with Q(v) <= u lie in
{Q <= (sqrt u + delta)^2}, delta = max sqrt Q(+-1/2, +-1/2), so N(u) := #{v != 0 : Q(v) <= u} <= A (sqrt u + delta)^2 - 1, A = 2 pi/sqrt|D|.
Partial summation: sum_{n>X} r(n) n^{-sigma} = -N(X) X^{-sigma} + sigma int_X^inf N(u) u^{-sigma-1} du
  <= -N(X) X^{-sigma} + A sigma [X^{1-sigma}/(sigma-1) + 2 delta X^{1/2-sigma}/(sigma-1/2) + delta^2 X^{-sigma}/sigma].
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import json, numpy as np
from math import sqrt
from flint import arb, ctx
ctx.prec = 80
FORMS=[(1,0,1),(1,0,5),(2,2,3),(1,1,6),(2,1,3),(1,0,14),(2,0,7),(3,2,5),(1,0,21),(3,0,7),(2,2,11),(5,4,5)]
X=10**6
def lam_of(Q):
    a,b,c=Q; return ((a+c)-sqrt((a-c)**2+b*b))/2
def reps(Q):
    a,b,c=Q; M=int(2*sqrt(X/lam_of(Q)))+3
    x,y=np.meshgrid(np.arange(-M,M+1),np.arange(-M,M+1))
    v=(a*x*x+b*x*y+c*y*y).ravel(); v=v[(v>0)&(v<=X)]
    return np.bincount(v, minlength=X+1)
def B_arb(Q,r,sig):  # rigorous (partial-sum lower bound, full upper bound incl. tail)
    ns=np.nonzero(r)[0]; m=int(ns[0]); rm=int(r[m]); s=arb(sig); S=arb(0)
    for n in ns[1:]: S+=arb(int(r[n]))*(arb(m)/arb(int(n)))**s
    S/=rm
    a,b,c=Q; D=b*b-4*a*c; A=2*arb.pi()/arb(-D).sqrt()
    dl=max(arb(a+b+c)/4, arb(a-b+c)/4).sqrt()
    NX=arb(int(r.sum())); Xa=arb(X)
    T=-NX*Xa**(-s)+A*s*(Xa**(1-s)/(s-1)+2*dl*Xa**(arb(0.5)-s)/(s-arb(0.5))+dl*dl*Xa**(-s)/s)
    tail=T*arb(m)**s/rm
    return S, S+tail
out=[]
for Q in FORMS:
    r=reps(Q); ns=np.nonzero(r)[0]; m=int(ns[0]); rm=int(r[m]); rr=r[ns[1:]]/rm; q=m/ns[1:]
    lo,hi=1.0,3.0
    for _ in range(50):    # float bisection; certification is the Arb step below
        mid=(lo+hi)/2
        if (rr*q**mid).sum()<1: hi=mid
        else: lo=mid
    dn=float(np.floor(lo*1e4)/1e4); ok_dn=B_arb(Q,r,dn)[0]>1
    up=float(np.ceil(hi*1e4)/1e4)
    for k in range(40):
        ok_up=B_arb(Q,r,up)[1]<1
        if ok_up: break
        up=round(up+0.0005,4)
    rec=dict(form=Q,m=m,r_m=rm,sigma_dom_lo=dn,sigma_dom_up=up,certified_up=bool(ok_up),certified_lo=bool(ok_dn))
    out.append(rec); print(rec,flush=True)
json.dump(out,open('dominant_abscissa.json','w'),indent=1)
