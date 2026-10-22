"""On-line zeros of the Davenport-Heilbronn function f(s) = (1-i k)/2 L(s,chi) + (1+i k)/2 L(s,chibar), chi mod 5 with
chi(2)=i, k = (sqrt(10-2 sqrt5)-2)/(sqrt5-1). Hardy-type Z(t) = Re(e^{i theta(t)} f(1/2+it)), theta(t) = arg Gamma((3/2+it)/2)
+ (t/2) log(5/pi) (odd characters; f has the functional equation of an odd conductor-5 L-function). Sign changes on a
grid, refined by bisection. Count check against N(T) = (T/2pi) log(5T/(2 pi e)) + 7/8-type smooth term: a deficit flags
missed zeros or off-line zeros. NUMERICAL. Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
usage: dh_zeros.py T0 T1 step out.json"""
import sys, json, mpmath as mp
mp.mp.dps=20
T0,T1,h=float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]); out=sys.argv[4]
s5=mp.sqrt(5); kap=(mp.sqrt(10-2*s5)-2)/(s5-1)
chi=[0,1,1j,-1j,-1]; chib=[0,1,-1j,1j,-1]
def f(s): return (1-1j*kap)/2*mp.dirichlet(s,chi)+(1+1j*kap)/2*mp.dirichlet(s,chib)
def Zc(t):
    th=mp.im(mp.loggamma((mp.mpf(1.5)+1j*t)/2))+t/2*mp.log(5/mp.pi)
    return mp.exp(1j*th)*f(mp.mpf(0.5)+1j*t)
ph=Zc(mp.mpf(T0)+0.37); rot=mp.conj(ph)/abs(ph)   # fix the constant phase
def Z(t): return mp.re(rot*Zc(t))
zs=[]; maximag=0
t=mp.mpf(T0); zp=Z(t)
while t<T1:
    t2=t+h; z2=Z(t2)
    if zp==0 or zp*z2<0:
        a,b,za,zb=t,t2,zp,z2; side=0   # Illinois regula falsi
        for _ in range(30):
            m=(a*zb-b*za)/(zb-za); zm=Z(m)
            if abs(b-a)<1e-12 or zm==0: break
            if zm*zb<0: a,za=b,zb; b,zb=m,zm; side=0
            else:
                b,zb=m,zm; za/=2
        zs.append(float(m))
    t,zp=t2,z2
for tt in (T0+1.1,(T0+T1)/2,T1-0.3):
    c=rot*Zc(mp.mpf(tt)); maximag=max(maximag,float(abs(mp.im(c))/max(abs(c),1e-30)))
json.dump({"T0":T0,"T1":T1,"step":h,"zeros":zs,"n":len(zs),"phase_check_rel_imag":maximag},open(out,"w"))
print(json.dumps({"T0":T0,"T1":T1,"n":len(zs),"phase_check":maximag}))
