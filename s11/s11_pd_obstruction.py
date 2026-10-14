# S11: witnesses |zeta(s+it)| > |zeta(s)| (violates PD bound |phi(t)|<=phi(0)) for sigma<=1
from mpmath import mp, zeta, mpf
mp.dps=30
for s in ['0.25','0.5','0.75','0.9','0.99','1']:
    s=mpf(s); z0=abs(zeta(s)) if s!=1 else None
    best=(0,0)
    for k in range(1,4000):
        t=mpf(k)/20; v=abs(zeta(s+1j*t))
        if v>best[1]: best=(t,v)
    print('sigma',s,'|zeta(sigma)|',z0,'max|zeta(sigma+it)| t<=200 at t=',best[0],best[1])
