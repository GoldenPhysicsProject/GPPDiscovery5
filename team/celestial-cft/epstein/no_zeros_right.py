"""CERTIFIED (exact rationals + rigorous tail bound): Z_Q(s) != 0 for Re s >= 3, for every form studied.
Let m = min Q over nonzero points and r(m) its multiplicity. Then
  |Z_Q(s) m^s / r(m) - 1| <= sum_{n>m} (r(n)/r(m)) (m/n)^sigma =: B(sigma),
which decreases in sigma. B(3) < 1 implies no zeros for sigma >= 3.
Tail: N(u) = #{v != 0 : Q(v) <= u} <= (2 sqrt(u/lam) + 1)^2 with lam = the min eigenvalue (box count), so
  sum_{n>X} r(n) n^{-sigma} <= sigma * int_X^inf N(u) u^{-sigma-1} du, bounded with N(u) <= C u for u >= X.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
from fractions import Fraction
from math import isqrt, sqrt
FORMS=[(1,0,1),(1,0,5),(2,2,3),(1,1,6),(2,1,3),(1,0,14),(2,0,7),(3,2,5),(1,0,21),(3,0,7),(2,2,11),(5,4,5)]
X=4000; sig=3
for Q in FORMS:
    a,b,c=Q; r={}
    M=isqrt(4*X)+3
    for x in range(-M,M+1):
        for y in range(-M,M+1):
            v=a*x*x+b*x*y+c*y*y
            if 0<v<=X: r[v]=r.get(v,0)+1
    m=min(r); rm=r[m]
    S=sum(Fraction(r[n],rm)*Fraction(m,n)**sig for n in r if n>m)
    lam=((a+c)-sqrt((a-c)**2+b*b))/2*0.999999
    C=(2/sqrt(lam)+1/sqrt(X))**2          # N(u) <= C u for u >= X  [(2 sqrt(u/lam)+1)^2 <= u (2/sqrt(lam)+1/sqrt(X))^2]
    tail=sig*C*X**(1-sig)/(sig-1)*m**sig/rm
    B=float(S)+tail
    print(Q,'m=',m,'r(m)=',rm,'B(3) <=',round(B,6),'OK' if B<1 else 'FAIL')
