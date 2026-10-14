# Hagedorn gap test: prime-gas partial sums Z_N(sigma)=sum_{n<=N} |a(n)| n^-sigma (float)
from mpmath import mp, mpf, log
mp.dps=20
N=200000
th=mpf('0.3')
def sig2(n,t):  # sum_{d|n} d^t * (n/d)^-t  = coeff of zeta(s-t)zeta(s+t)
    s=0;d=1
    while d*d<=n:
        if n%d==0:
            e=n//d; s+= (d/e)**t+(e/d)**t if d!=e else 1
        d+=1
    return s
for sig in [mpf('1.1'),mpf('1.2'),mpf('1.3')]:
    for M in [10**3,10**4,10**5]:
        z=sum(mpf(n)**-sig for n in range(1,M+1))
        f=sum(sig2(n,float(th))*mpf(n)**-sig for n in range(1,M+1)) if M<=10**4 else None
        print(f"sigma={sig} N={M} Z_zeta={float(z):.6f} Z_Ftheta={f if f is None else round(float(f),3)}")
