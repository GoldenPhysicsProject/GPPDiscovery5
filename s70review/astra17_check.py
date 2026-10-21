"""Grok S70 review: independent (NON-certified, mpmath) evaluation of the zeta Weil form on Astra's
17-mode Fourier block e^{i pi j x/L}/sqrt(2L), j=-8..8, on [-L,L], using Grok's S14/S62 Weil functional
(same formula as s62/cosritz.py) but a fresh quadrature implementation. Split by parity:
even span = cos(j pi x/L) j=0..8 ; odd span = sin(j pi x/L) j=1..8.  Usage: astra17_check.py L dps"""
import sys, mpmath as mp
L=mp.mpf(sys.argv[1]); mp.mp.dps=int(sys.argv[2]); J=8
def mang(n):
    for p in range(2,n+1):
        if n%p==0:
            m=n
            while m%p==0: m//=p
            return p if m==1 else None
ns=[n for n in range(2,int(mp.exp(2*L))+2) if mang(n) and mp.log(n)<2*L]
c0=-mp.euler-mp.log(mp.pi)-mp.log(1-mp.exp(-4*L))
k1=lambda x: mp.exp(-x)/mp.sinh(x); k2=lambda x: 2*mp.exp(-x/4)*mp.sinh(3*x/4)/mp.sinh(x)
def icos(k,c,al,be):  # int_al^be cos(k y + c) dy
    return (be-al)*mp.cos(c) if k==0 else (mp.sin(k*be+c)-mp.sin(k*al+c))/k
def F(par,am,an,x):   # int_{x-L}^{L} f_m(y) f_n(y-x) dy, x in [0,2L]
    P=icos(am+an,-an*x,x-L,L); M=icos(am-an,an*x,x-L,L)
    return (P+M)/2 if par=='even' else (M-P)/2
def W(par,am,an):
    Fs=lambda x:(F(par,am,an,x)+F(par,an,am,x))/2
    F0=Fs(mp.mpf(0)); seg=[0,L/4,L/2,L,3*L/2,2*L]
    r=c0*F0+mp.quad(lambda x:(F0-Fs(x))*k1(x),seg)-mp.quad(lambda x:Fs(x)*k2(x),seg)
    r-=sum(2*mp.log(mang(n))/mp.sqrt(n)*Fs(mp.log(n)) for n in ns)
    return r
def pole(par,a):
    if par=='even': return mp.quad(lambda u:mp.cos(a*u)*mp.cosh(u/2),[-L,0,L])
    return mp.quad(lambda u:mp.sin(a*u)*mp.sinh(u/2),[-L,0,L])
for par in ['even','odd']:
    js=list(range(0,J+1)) if par=='even' else list(range(1,J+1))
    a=[j*mp.pi/L for j in js]; K=len(a); sg=2 if par=='even' else -2
    p=[pole(par,x) for x in a]
    Wm=mp.matrix(K,K)
    for i in range(K):
        for k in range(i,K):
            Wm[i,k]=Wm[k,i]=W(par,a[i],a[k])+sg*p[i]*p[k]
    g=[mp.sqrt(2*L if (par=='even' and js[i]==0) else L) for i in range(K)]
    N=mp.matrix(K,K)
    for i in range(K):
        for k in range(K): N[i,k]=Wm[i,k]/(g[i]*g[k])
    E=mp.eigsy(N,eigvals_only=True)
    print(par,"L=",sys.argv[1],"min=",mp.nstr(min(E),15),"next=",mp.nstr(sorted(E)[1],8),flush=True)
