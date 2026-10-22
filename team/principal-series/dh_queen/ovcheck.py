import json, mpmath as mp, sys
from fractions import Fraction
mp.mp.dps=40
g0,b0=mp.mpf('85.699348'),mp.mpf('0.308517')
fn=sys.argv[1]; d=json.load(open(fn)); L=mp.mpf(Fraction(d['L'].replace('_','/')).numerator)/Fraction(d['L'].replace('_','/')).denominator
v=[mp.mpf(Fraction(s).numerator)/Fraction(s).denominator for s in d['vec']]
# exact: int_{-L}^{L} cos(a x) sin(g x) sinh(b x) dx = 0 (odd integrand!) -- so use the EVEN Q part: -sin(g x)sinh(b x) is even; cos*even -> even
a=[(2*m+1)*mp.pi/(2*L) for m in range(len(v))]
def I(am):  # 2 int_0^L cos(am x) sin(g0 x) sinh(b0 x) dx, closed form via exponentials
    tot=mp.mpc(0)
    for s1 in (1,-1):
        for s2 in (1,-1):
            for s3 in (1,-1):
                # cos= (e^{i am x}+e^{-i am x})/2, sin=(e^{igx}-e^{-igx})/(2i), sinh=(e^{bx}-e^{-bx})/2
                c=mp.mpf(1)/8/mp.mpc(0,1)*(s2)*(s3)
                k=mp.mpc(s3*b0, s1*am+s2*g0)
                tot+=c*(mp.e**(k*L)-1)/k
    return 2*mp.re(tot)
fg=mp.fsum(vi*I(ai) for vi,ai in zip(v,a))
nf=mp.sqrt(L*mp.fsum(t*t for t in v))
ng=mp.sqrt(mp.quad(lambda x:(mp.sin(g0*x)*mp.sinh(b0*x))**2,mp.linspace(-L,L,200)))
fpL=mp.fsum(-vi*ai*mp.sin(ai*L) for vi,ai in zip(v,a))
print(d['L'],len(v),'overlap(Im Q)=',mp.nstr(fg/nf/ng,8),"f'(L)/|f|=",mp.nstr(fpL/nf,8),'f(0)/|f|=',mp.nstr(mp.fsum(v)/nf,8))
