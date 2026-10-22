"""[Goddess copy of GPPDiscovery5/s62/cosritz.py (Grok S62), assembly unchanged; tail replaced by a matrix dump]
Grok S62 H62a: rigorous L2 Ritz UPPER bounds for the zeta Weil form on even (and odd) L2[-L,L],
Zhu-type trig basis (same span as Zhu arXiv:2608.24827 Sec.11 sine modes), geometric side only, Arb.
EVEN modes c_m(u)=cos(a_m u), a_m=(2m+1)pi/(2L); ODD modes s_m(u)=sin(b_m u), b_m=m pi/L (m>=1). Gram = L*I.
Autocorrelation F_mn(x)=int c_m(y)c_n(y-x)dy on [0,2L] (derived in closed form, eps=(-1)^(m+n)):
  even m!=n: eps[a_m sin(a_n x) - a_n sin(a_m x)]/(a_m^2-a_n^2);  m=n: (2L-x)cos(a_m x)/2 + sin(a_m x)/(2a_m)
  odd s_m: identical closed forms with a->b (verified numerically, selftest.py)
Weil functional on even F (S14 normalisation, validated S27/S39/S61):
  W(F) = pole + c0 F(0) + int_0^{2L}(F0-F)k1 - int_0^{2L} F k2 - sum_{log n<2L} 2Lambda(n)/sqrt(n) F(log n),
  k1=e^{-x}/sinh x, k2=2e^{-x/4}sinh(3x/4)/sinh x, c0=-gamma-log pi-log(1-e^{-4L});
  pole = 2(int f cosh(u/2))^2 (even) / -2(int f sinh(u/2))^2 (odd).
All integrals: acb.integral (rigorous), integrands rewritten analytic (sinc). Ritz vector (mpmath) -> rationals ->
Arb Rayleigh quotient v^T W v / (L |v|^2) = rigorous upper bound on inf Q/||f||^2 in the sector (min-max).
Usage: cosritz.py L K prec dps [even|odd|both]"""
import sys, json, time
from flint import arb, acb, arb_mat, fmpq, ctx
import mpmath as mp
Ls=sys.argv[1]; K=int(sys.argv[2]); prec=int(sys.argv[3]); dps=int(sys.argv[4]); sec=sys.argv[5] if len(sys.argv)>5 else "even"
ctx.prec=prec; mp.mp.dps=dps
def parseL(s):
    if '/' in s: p,q=s.split('/'); return arb(int(p))/int(q)
    return arb(s)
L=parseL(Ls); I=acb(0,1); PI=arb.pi()
def mang(n):
    for p in range(2,n+1):
        if n%p==0:
            m=n
            while m%p==0: m//=p
            return p if m==1 else None
ns=[n for n in range(2,int(float((2*L).exp().mid()))+2) if mang(n) and arb(n).log()<2*L]
assert all(not (arb(n).log()-2*L).contains(0) for n in ns)
w={n:2*arb(mang(n)).log()/arb(n).sqrt() for n in ns}; lg={n:arb(n).log() for n in ns}
c0=-arb.const_euler()-PI.log()-(1-(-4*L).exp()).log()
import os
FORM=os.environ.get('FORM','ZETA')   # FORM=DH: Davenport-Heilbronn form (S47 difference trick; no pole, +log5 F(0) + int_0^{2L} F sech(x/2), b(n) comb)
if FORM=='DH':
    def chi5(n): return {0:(0,0),1:(1,0),2:(0,1),3:(0,-1),4:(-1,0)}[n%5]
    s5=arb(5).sqrt(); kap=((10-2*s5).sqrt()-2)/(s5-1)
    NN=int(float((2*L).exp().mid())); nn=[n for n in range(2,NN+2) if arb(n).log()<2*L]
    aa={n:arb(chi5(n)[0])+kap*chi5(n)[1] for n in range(1,max(nn)+1)}; bb={1:arb(0)}
    for n in range(2,max(nn)+1): bb[n]=aa[n]*arb(n).log()-sum((bb[d]*aa[n//d] for d in range(2,n) if n%d==0),arb(0))
    ns=nn; w={n:2*bb[n]/arb(n).sqrt() for n in ns}; lg={n:arb(n).log() for n in ns}
def integ(f): return acb.integral(lambda z,an: f(z),0,2*L).real
def k1(z): return (-z).exp()/(I*z).sinc()            # times 1/x : k1 = e^{-x}/sinh x = that/x
def k2(z): return 2*(-z/4).exp()*(I*z*3/4).sinc()*arb(3)/4/(I*z).sinc()
def Lsin(a):   # W'(sin a x) without pole
    A3=integ(lambda z: a*(a*z).sinc()*k1(z)); B3=integ(lambda z:(a*z).sin()*k2(z))
    r=-A3-B3-sum(w[n]*(a*lg[n]).sin() for n in ns)
    if FORM=='DH': r+=integ(lambda z:(a*z).sin()/(z/2).cosh())
    return r
def LCx(a):    # W'((2L-x)cos(a x)/2)
    A1=integ(lambda z: a*(a*z/2).sin()*(a*z/2).sinc()*k1(z)); A2=integ(lambda z:(a*z).cos()*k1(z))
    B1=integ(lambda z:(a*z).cos()*k2(z)); B2=integ(lambda z:z*(a*z).cos()*k2(z))
    r=c0*L+L*A1+A2/2-L*B1+B2/2-sum(w[n]*(2*L-lg[n])*(a*lg[n]).cos()/2 for n in ns)
    if FORM=='DH': r+=arb(5).log()*L+integ(lambda z:(2*L-z)*(a*z).cos()/2/(z/2).cosh())
    return r
t0=time.time(); out={"form":FORM,"L":Ls,"K":K,"prec":prec,"dps":dps,"ns":ns}
def assemble_even():
    a=[(2*m+1)*PI/(2*L) for m in range(K)]
    S=[Lsin(x) for x in a]; C=[LCx(x) for x in a]
    PF=0 if FORM=='DH' else 1
    cp=[PF*2*a[m]*(-1)**m*(L/2).cosh()/(a[m]**2+arb(1)/4) for m in range(K)]
    W=arb_mat(K,K)
    for m in range(K):
        for n in range(K):
            if m==n: W[m,m]=2*cp[m]**2+C[m]+S[m]/(2*a[m])
            else: W[m,n]=2*cp[m]*cp[n]+(-1)**(m+n)*(a[m]*S[n]-a[n]*S[m])/(a[m]**2-a[n]**2)
    return W
def assemble_odd():
    # s_m(u)=sin(b_m u), b_m=m pi/L, m=1..K. F_mn(x)=int_{x-L}^{L} sin(b_m y) sin(b_n(y-x)) dy (x in [0,2L]):
    # m!=n: eps[b_m sin(b_n x)-b_n sin(b_m x)]/(b_m^2-b_n^2); m=n: (2L-x)cos(b x)/2 + sin(b x)/(2b)  (same form as even; selftest.py)
    b=[m*PI/L for m in range(1,K+1)]
    S=[Lsin(x) for x in b]; C=[LCx(x) for x in b]
    sp=[]   # int_{-L}^{L} sin(b u) sinh(u/2) du = 2 int_0^L ... = 2[(-b cos(bu)sinh(u/2)+ sin(bu)cosh(u/2)/2)/(b^2+1/4)]_0^L
    for m in range(K): sp.append((0 if FORM=='DH' else 1)*2*(-b[m]*(-1)**(m+1)*(L/2).sinh())/(b[m]**2+arb(1)/4))
    W=arb_mat(K,K)
    for m in range(K):
        for n in range(K):
            if m==n: W[m,m]=-2*sp[m]**2+C[m]+S[m]/(2*b[m])
            else: W[m,n]=-2*sp[m]*sp[n]+(-1)**(m+n)*(b[m]*S[n]-b[n]*S[m])/(b[m]**2-b[n]**2)
    return W

def dump(W,nm):
    rows=[[ (W[i,j].mid().str(radius=False, more=True), float(W[i,j].rad())) for j in range(K)] for i in range(K)]
    json.dump({"form":FORM,"L":Ls,"K":K,"sector":nm,"prec":prec,"ns":ns,"W":rows,"sec":round(time.time()-t0,1)},open(f"Wdump_{FORM}_{nm}_L{Ls.replace('/','_')}_K{K}_p{prec}.json","w"))
if sec in("even","both"): dump(assemble_even(),"even")
if sec in("odd","both"): dump(assemble_odd(),"odd")
print(json.dumps({"done":True,"form":FORM,"K":K,"sec":round(time.time()-t0,1)}),flush=True)
