"""NUMERICAL hunt for zeros of Z_{x^2+5y^2} with Re s > 1 (Kronecker / Bohr screening). The result is then CERTIFIED by
certify_gt1.py. For Re s > 1, Z = L4 L5 (1 + R), with
R(s) = (1+2^-s)/(1-2^-s) * prod_{p = 3,7 mod 20} ((1+p^-s)/(1-p^-s))^2   (zeta L_-20 / (L_-4 L_5), from the Euler factors).
Stage 1 keeps t whose 2-, 3-, 7-adic phases are favourable. Stage 2 uses the truncated product over p < 5000 and keeps local minima of |1+R|.
Stage 3 runs Newton on the exact 1 + R (Arb Hurwitz L-values).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, numpy as np, math
sys.path.insert(0,'.')
from sympy import primerange
from zq_L import Ld
from flint import acb, ctx
ctx.prec=64
SIG=float(sys.argv[1]); T0=float(sys.argv[2]); T1=float(sys.argv[3]); h=0.004
ps=np.array([p for p in primerange(3,5000) if p%20 in (3,7)],float)
def logR(s, P):
    z2=2.0**(-s); out=np.log((1+z2)/(1-z2))
    for p in P:
        z=p**(-s); out=out+2*np.log((1+z)/(1-z))
    return out
def Gexact(s):
    s=acb(s.real,s.imag); v=1+s.zeta()*Ld(-20,s)/(Ld(-4,s)*Ld(5,s))
    return complex(float(v.real.mid()),float(v.imag.mid()))
def newton(s):
    for _ in range(40):
        g=Gexact(s); d=(Gexact(s+1e-7)-g)/1e-7
        if d==0: return None
        step=g/d; s=s-step
        if abs(step)<1e-13: return s
        if not (0.9<s.real<2) : return None
    return None
found=[]
for a in np.arange(T0,T1,1e4):
    t=np.arange(a,min(a+1e4,T1),h); s=SIG+1j*t
    l3=logR(s,ps[:2])                        # p = 2, 3, 7
    m=(np.abs(l3.imag)>2.0)&(np.abs(l3.real)<0.8)
    if not m.any(): continue
    tc=t[m]; L=logR(SIG+1j*tc,ps)
    v=np.abs(1+np.exp(L))
    # local minima among contiguous candidates
    for i in range(1,len(tc)-1):
        if v[i]<0.15 and v[i]<=v[i-1] and v[i]<=v[i+1]:
            z=newton(complex(SIG,tc[i]))
            if z is not None and z.real>1 and all(abs(z-f)>1e-6 for f in found):
                found.append(z); print(json.dumps(dict(x=z.real,y=z.imag,seed_t=float(tc[i]),trunc_abs=float(v[i]))),flush=True)
print('done',T0,T1,flush=True)
