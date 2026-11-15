"""(T2-rig): exact reduction of omega^4 rho_j to Q(omega)=R(omega)+Pi(omega); numerical test against Round 8 values."""
import numpy as np, json
from scipy import special
l=0.25; NS=[(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]; bI2=16*np.sinh(l)**2
om=lambda j: j*np.pi/l
sig=lambda j: -np.sqrt(2/l)*(-1)**j*om(j)*np.sinh(l/2)/(0.25+om(j)**2)
# P(t)=sin^4(t/2) sin^2(t/4) = sum_k c_k cos(k t/2)
M=4096; tt=np.arange(M)*4*np.pi/M; Pv=np.sin(tt/2)**4*np.sin(tt/4)**2
F=np.fft.rfft(Pv)/M; c=np.real(F[:6]); c[1:]*=2
c=np.round(c*1024)/1024  # exact dyadic
assert np.max(np.abs(sum(c[k]*np.cos(k*tt/2) for k in range(6))-Pv))<1e-12
Cn={n:float(sum(c[k] for k in range(6) if k/2<np.log(n))) for n,_ in NS}
A=lambda t: special.digamma(0.25+0.5j*t).real
x0,w0=np.polynomial.legendre.leggauss(40)
def R(w,T=4e5,h=1.0):
    # integrand smooth (removable at t=w); panels
    ed=np.arange(0,T,h); t=(ed[:,None]+h/2*(x0[None,:]+1)).ravel(); wt=np.tile(h/2*w0,len(ed))
    P=np.sin(t/2)**4*np.sin(t/4)**2
    d=t*t-w*w; g=np.where(np.abs(d)>1e-9,(A(t)-A(w))/np.where(np.abs(d)>1e-9,d,1),0)
    val=w*np.sum(wt*g*P)
    tail=w*(3/16)*(np.log(T/w)+1)/T   # mean-P tail of log(t/w)/t^2
    return val+tail
def Pi(w): return np.pi*sum(p_ := 0 for _ in [])+np.pi*sum(np.log(p)/np.sqrt(n)*Cn[n]*np.sin(w*np.log(n)) for n,p in NS)
kap=9/(4*np.sqrt(3)); A0=A(0.0)
def Rbar(w): return np.pi**2/4+kap/w*np.log1p(w)+w*(np.log(w)-np.log(2)+kap/w-A0)/(w*w-1)
Pib=np.pi*sum(np.log(p)/np.sqrt(n)*abs(Cn[n]) for n,p in NS)
pref=(16/np.pi)*np.sqrt(2/l)
out=dict(c=c.tolist(),Cn=Cn,kappa=kap,A0=float(A0),Pi_bound=float(Pib),checks=[])
ref={73:-5.35066384039753,123:-3.101766731280355,173:-2.8002357308944315,273:-12.956305194499834,472:8.055786994250864}
for j in [13,73,123,173,273,472,1012]:
    w=om(j); Rv=R(w); Pv_=Pi(w); Q=Rv+Pv_
    val=pref*(-1)**j*Q-(bI2**2/8)*np.sinh(l/2)*sig(j)
    out['checks'].append(dict(j=j,R=float(Rv),Pi=float(Pv_),Q=float(Q),w4rho_delta=float(val),ref=ref.get(j),Rbar=float(Rbar(w))))
    print(j,Rv,Pv_,val,ref.get(j),Rbar(w),flush=True)
print('c',c,'Cn',Cn,'Pi_bound',Pib,'Rbar(om13)',Rbar(om(13)),'Rbar(4pi)',Rbar(om(1)))
json.dump(out,open('t2rig_results.json','w'),indent=1)
