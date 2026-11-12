# Turbo 6821018808 item 5: H for zeta, triangle psi. Zero-free (prime/arch) side vs zero side of K(d).
import numpy as np, mpmath as mp
from scipy.special import digamma
from scipy.integrate import quad
from sympy import factorint
import json, sys
def Psi(t):
    t=np.asarray(t,float); x=t/2; s=np.where(np.abs(x)<1e-8,1.0,np.sin(x)/np.where(x==0,1,x)); return s**2
def g(y):  # psi*psi, psi = max(0,1-|x|)
    y=abs(y)
    if y<=1: return 2/3 - y*y + y**3/2
    if y<=2: return (2-y)**3/6
    return 0.0
def Lam(n):
    f=factorint(n)
    return float(np.log(list(f)[0])) if len(f)==1 else 0.0
NMAX=int(np.exp(8))+2
LAM=[0,0]+[Lam(n) for n in range(2,NMAX)]
Psi_i2=(float(mp.sinh(0.25))/0.25)**2   # Psi(i/2)=Psi(-i/2)
def arch(d):
    f=lambda t: Psi(t)*Psi(t)*np.cos(t*d)*(digamma(0.25+0.5j*t).real-np.log(np.pi))/(2*np.pi)
    tot=0.0; a=0.0
    for b in list(np.arange(20,4000,20)):
        tot+=quad(f,a,b,limit=400,epsabs=1e-14)[0]; a=b
    return 2*tot   # integrand even in t
def K_prime(d):
    pole=2*Psi_i2**2*np.cosh(d/2)
    P=sum(LAM[n]/np.sqrt(n)*(g(np.log(n)+d)+g(-np.log(n)+d)) for n in range(2,NMAX))
    return pole-P+arch(d)
mp.mp.dps=20
NZ=int(sys.argv[1]) if len(sys.argv)>1 else 400
gam=np.array([float(mp.zetazero(k).imag) for k in range(1,NZ+1)])
T=gam[-1]
tail=256*(np.log(T/(2*np.pi))+1)/(2*np.pi*3*T**3)  # bound on sum_{gamma>T} 2Psi^2 (|cos|<=1), crude
def K_zero(d): return float(np.sum(2*Psi(gam)**2*np.cos(gam*d)))
out={"NZ":NZ,"T":T,"tail_bound":2*tail,"Psi_i2":Psi_i2,"rows":[]}
for d in [0,0.25,0.5,1,1.5,2,3,4,5,6]:
    kp=K_prime(d); kz=K_zero(d)
    out["rows"].append({"d":d,"K_prime_side":kp,"K_zero_side":kz,"diff":kp-kz})
    print(d,kp,kz,kp-kz,flush=True)
json.dump(out,open("hef_tri.json","w"),indent=1)
