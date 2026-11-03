# NUMERICAL (float64) scan for P.17.0(i): smallest T at which (H3^sigma)(T, eta, mu) admits mu > 0 at 2L = log 3. Not a certificate. RH OPEN.
import numpy as np
from scipy.linalg import eigh
exec(open('p17_h3_probe.py').read().split("for par in")[0])
c3=b/np.sqrt(3); c0=(c3+np.sqrt(c3**2+4*c2**2))/2
def mu_max(H0,Nlo,G):
    lo,hi_=-3.0,2.0
    if eigh(H0-lo*Nlo,G,eigvals_only=True)[0]<0: return None
    for _ in range(50):
        mid=(lo+hi_)/2
        if eigh(H0-mid*Nlo,G,eigvals_only=True)[0]>=0: lo=mid
        else: hi_=mid
    return lo
N=40
for par in ['odd','even']:
    ns=[n for n in range(1,4*N) if (n%2==0)==(par=='odd')][:N]
    t=np.linspace(0,4000,400001)[1:]; dt=t[1]-t[0]
    F=np.array([basis_ft(n,t) for n in ns]); w=(0.25+t*t)/np.pi*dt
    gram=lambda s: np.real((F*(w*s))@F.conj().T)
    G=gram(np.ones_like(t)); Q2=gram(Om(t)-2*c2*np.cos(a*t))
    best=[]
    for T in [16,20,25,30,35,38.1,45,55,70]:
        Nhi=gram((t>=T)*1.0); Nlo=G-Nhi
        row=[]
        for eta in [0.0373,0.1,0.25,0.5,1,2,4]:
            if c0*(1+eta)>Om(T)+1e-12: row.append(None); continue
            m=mu_max(Q2-c0*(1+eta)*Nhi,Nlo,G)
            l3=None if (m is None or m<=0) else np.pi*m/(4*c0*(1+1/eta)*T)
            row.append((eta,m,l3))
        print(par,'T=%5.1f Om=%.3f'%(T,Om(T)),' '.join('eta=%g:mu=%s,l3=%s'%(r[0],'%.4g'%r[1] if r[1] is not None else 'fail','%.2e'%r[2] if r[2] else '-') for r in row if r))
