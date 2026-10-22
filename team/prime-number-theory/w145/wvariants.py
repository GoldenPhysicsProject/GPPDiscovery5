import numpy as np, math
from mpmath import euler
g=float(euler); N=2*10**6; T=14.5
lam=np.zeros(N+1); sv=np.ones(N+1,bool); sv[:2]=False
for p in range(2,int(N**.5)+1):
    if sv[p]: sv[p*p::p]=False
pr=np.nonzero(sv)[0]
lamP=np.zeros(N+1); lamP[pr]=np.log(pr)
for p in pr:
    q=p
    while q<=N: lam[q]=math.log(p); q*=p
def grid_W(L,h,rule):
    # Astra §3 on a uniform t-grid: j-c = -2 int_t^Tmax q - int_0^t e^{-(t-u)/2} q ; r=e^{t/2}(j-c); W=int r^2
    psi=np.cumsum(L)
    Tm=math.log(N); t=np.arange(0,Tm+h/2,h)
    if rule=='left': us=t[:-1]
    else: us=t[:-1]+h/2
    qv=psi[np.floor(np.exp(us)).astype(int)]*np.exp(-us)-1
    # future: int_t^Tm q
    cum=np.concatenate(([0],np.cumsum(qv*h))); fut=cum[-1]-cum
    # past damped: recursion
    past=np.zeros_like(t)
    for k in range(1,len(t)):
        past[k]=past[k-1]*math.exp(-h/2)+qv[k-1]*h*math.exp(-h/4)
    r=np.exp(t/2)*(-2*fut-past)
    m=t<=T
    return np.trapezoid(r[m]**2,t[m])
for h in [1e-2,1e-3]:
    for rule in ['left','mid']:
        print("psi grid h=%g %s: W(14.5)=%.3f"%(h,rule,grid_W(lam,h,rule)))
print("theta (primes only) grid h=1e-3 mid:",grid_W(lamP,1e-3,'mid'))
print('---scan')
for h in [0.005,0.01,0.02,0.025,0.05,0.1]:
    print(h, round(grid_W(lam,h,'mid'),2), round(grid_W(lam,h,'left'),2))
