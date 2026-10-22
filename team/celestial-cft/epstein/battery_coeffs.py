"""Battery hook: Dirichlet coefficients c(n) of -Z_Q'/Z_Q (principal forms, a=1, so Z_Q/r(1) = 1 + sum_{n>=2} r(n)/r(1) n^{-s})
compared with those of the Dedekind zeta zeta_K, K = Q(sqrt D), which has the same Gamma factor Gamma(s)(sqrt|D|/2pi)^s.
For zeta_K, c(n) = Lambda_K(n) >= 0 and is supported on prime powers. For Z_Q, with no Euler product, we look at signs and support.
Exact rational arithmetic.  Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json
from fractions import Fraction
from math import log, isqrt
from sympy import factorint, kronecker_symbol

def r_Q(Q,N):
    a,b,c=Q; r=[0]*(N+1); M=isqrt(4*N)+3
    for x in range(-M,M+1):
        for y in range(-M,M+1):
            v=a*x*x+b*x*y+c*y*y
            if 0<v<=N: r[v]+=1
    return r

def logderiv_coeffs(f,N):
    # f[1]=1; g = -f'/f has coefficients c(n) with sum_{d|n} c(d) f(n/d) = f(n) log n.  Returns c(n)/log-free form:
    # c(n) = f(n) log n - sum_{d|n, d<n} c(d) f(n/d)   (floats)
    c=[0.0]*(N+1)
    for n in range(2,N+1):
        s=f[n]*log(n)
        for d in range(2,n):
            if n%d==0: s-=c[d]*f[n//d]
        c[n]=s
    return c

N=int(sys.argv[1]) if len(sys.argv)>1 else 200
out={}
for Q in [(1,0,5),(1,1,6),(1,0,14),(1,0,21)]:
    a,b,cc=Q; D=b*b-4*a*cc
    r=r_Q(Q,N); f=[0.0]+[r[n]/r[1] for n in range(1,N+1)]
    cQ=logderiv_coeffs(f,N)
    # Dedekind zeta_K = zeta * L(chi_D): coefficients sum_{d|n} chi_D(d)
    fK=[0.0]+[sum(int(kronecker_symbol(D,d)) for d in range(1,n+1) if n%d==0) for n in range(1,N+1)]
    cK=logderiv_coeffs(fK,N)
    nonpp=[n for n in range(2,N+1) if abs(cQ[n])>1e-9 and len(factorint(n))>1]
    neg=[n for n in range(2,N+1) if cQ[n]<-1e-9]
    out[str(Q)]=dict(D=D,first_negative=neg[:10],n_negative=len(neg),first_non_prime_power_support=nonpp[:10],
        sample=[(n,round(cQ[n],6),round(cK[n],6)) for n in range(2,31)])
    print(Q,'D=',D,'#neg c(n), n<=%d:'%N,len(neg),'first neg',neg[:8],'non-prime-power support',nonpp[:8])
json.dump(out,open('battery_coeffs.json','w'),indent=1)
