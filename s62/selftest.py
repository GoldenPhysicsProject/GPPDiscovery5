# numeric check of the closed-form autocorrelations F_mn(x)=int_{x-L}^{L} f_m(y) f_n(y-x) dy, and of the pole moments
import mpmath as mp; mp.mp.dps=30
L=mp.mpf('1.35'); err=0
for kind in ("even","odd"):
    fr=(lambda m:(2*m+1)*mp.pi/(2*L)) if kind=="even" else (lambda m:(m+1)*mp.pi/L)
    f=(lambda a,y: mp.cos(a*y)) if kind=="even" else (lambda a,y: mp.sin(a*y))
    for (m,n) in [(0,0),(3,3),(0,1),(2,5),(4,1)]:
        a,b=fr(m),fr(n); eps=(-1)**(m+n)
        for x in [mp.mpf('0.3'),mp.mpf('1.1'),mp.mpf('2.2')]:
            num=mp.quad(lambda y: f(a,y)*f(b,y-x),[x-L,L])
            cf=((2*L-x)*mp.cos(a*x)/2+mp.sin(a*x)/(2*a)) if m==n else eps*(a*mp.sin(b*x)-b*mp.sin(a*x))/(a**2-b**2)
            err=max(err,abs(num-cf))
    for m in range(4):
        a=fr(m)
        if kind=="even": num=mp.quad(lambda u: mp.cos(a*u)*mp.cosh(u/2),[-L,L]); cf=2*a*(-1)**m*mp.cosh(L/2)/(a**2+mp.mpf(1)/4)
        else: num=mp.quad(lambda u: mp.sin(a*u)*mp.sinh(u/2),[-L,L]); cf=2*(-a*(-1)**(m+1)*mp.sinh(L/2))/(a**2+mp.mpf(1)/4)
        err=max(err,abs(num-cf))
print("max abs err closed forms:",mp.nstr(err,5))
