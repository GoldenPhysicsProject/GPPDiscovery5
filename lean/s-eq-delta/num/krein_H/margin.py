# Turbo 6821237527 item 2: margin decay of min Q on ker l_c (triangle psi, zero-free side), extended grid.
import numpy as np, json
from scipy.special import digamma
from scipy.integrate import quad
import sys
h=0.25; KMAX=int(sys.argv[1]) if len(sys.argv)>1 else 40
DMAX=KMAX*h; NMAX=int(np.exp(DMAX+2))+2
# von Mangoldt sieve
LAM=np.zeros(NMAX)
sieve=np.ones(NMAX,bool); sieve[:2]=False
for p in range(2,int(NMAX**.5)+1):
    if sieve[p]: sieve[p*p::p]=False
for p in np.nonzero(sieve)[0]:
    q=p
    while q<NMAX: LAM[q]=np.log(p); q*=p
ns=np.arange(NMAX); w=np.zeros(NMAX); w[2:]=LAM[2:]/np.sqrt(ns[2:]); ln=np.zeros(NMAX); ln[2:]=np.log(ns[2:])
def g(y):
    y=np.abs(y); return np.where(y<=1,2/3-y*y+y**3/2,np.where(y<=2,(2-y)**3/6,0.0))
def Psi(t):
    x=t/2; return 1.0 if abs(x)<1e-8 else (np.sin(x)/x)**2
Psi_i2=(np.sinh(0.25)/0.25)**2
def arch(d):
    f=lambda t: Psi(t)**2*np.cos(t*d)*(digamma(0.25+0.5j*t).real-np.log(np.pi))/(2*np.pi)
    tot=0.0;a=0.0
    for b in np.arange(20,4000,20): tot+=quad(f,a,b,limit=400,epsabs=1e-14)[0]; a=b
    return 2*tot
def K(d): return 2*Psi_i2**2*np.cosh(d/2)-np.sum(w*(g(ln+d)+g(-ln+d)))+arch(d)
Kp=np.array([K(k*h) for k in range(KMAX+1)])
old=np.array(json.load(open("scan.json"))["Kp"]); print("max |Kp - old Kp| on k<=24:",np.max(np.abs(Kp[:25]-old)))
def run(n,s):
    x=(np.arange(n)-(n-1)/2)*h+s
    M=Kp[np.abs(np.subtract.outer(np.arange(n),np.arange(n)))]
    a=np.cosh(x/2);a/=np.linalg.norm(a)
    U,_,_=np.linalg.svd(np.eye(n)-np.outer(a,a));B=U[:,:n-1]
    return float(np.linalg.eigvalsh(B.T@M@B)[0]), float(np.linalg.eigvalsh(M)[0])
rows=[];best=[]
for n in range(2,KMAX+2):
    bm=None
    for s in np.arange(0,3.01,0.25):
        m,mf=run(n,s); L=1+(n-1)*h/2+s
        rows.append({"n":n,"shift":float(s),"L":float(L),"min_Q_on_ker_lc":m,"min_eig_full":mf})
        if bm is None or m<bm["min_Q_on_ker_lc"]: bm=rows[-1]
    best.append(bm); print(n,bm["shift"],bm["L"],bm["min_Q_on_ker_lc"],bm["min_eig_full"],flush=True)
json.dump({"h":h,"Kp":Kp.tolist(),"rows":rows,"best_per_n":best},open("margin.json","w"),indent=1)
