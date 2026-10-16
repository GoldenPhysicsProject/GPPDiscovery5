"""Grok S42 H42b: audit Astra gamma-free Poisson transport (PR11 6215916305) + Grok hypothesis:
the Legendre-optimal (common-tilt, fully correlated) coupling dominates EVERY coupling on Astra's reference space,
so the coupling problem reduces exactly to the piecewise-constant inequality C_k >= V*(B_k). Tested to x=1e7 (float64 => NUMERICAL).
RH OPEN."""
import sympy as sp, mpmath as mp, numpy as np, json, math
out={}
t,lam,a,u=sp.symbols('t lambda alpha u',positive=True)
# 1. cgf of (N+ - N-)/2, N+-~Poisson(4): log E e^{tY} = 4(e^{t/2}-1)+4(e^{-t/2}-1)
mp.mp.dps=40
def cgf_num(tv,lamv=4,M=200):
    s=mp.nsum(lambda k: mp.e**(-lamv)*lamv**k/mp.factorial(k)*mp.e**(tv*k/2),[0,mp.inf])
    s2=mp.nsum(lambda k: mp.e**(-lamv)*lamv**k/mp.factorial(k)*mp.e**(-tv*k/2),[0,mp.inf])
    return mp.log(s*s2)
out['cgf_check']=[(tv, mp.nstr(cgf_num(mp.mpf(tv))-8*(mp.cosh(mp.mpf(tv)/2)-1),5)) for tv in ['0.3','1','2.5','7']]
# 2. alpha_n = w_n / sinh(u_n/2) = 2 Lambda(n)/(n-1)
n=sp.symbols('n',positive=True)
out['alpha_identity']=sp.simplify((sp.log(n)/sp.sqrt(n)/sp.sinh(sp.log(n)/2)).rewrite(sp.exp) - 2*sp.log(n)/(n-1))==0
# 3. tilt mean and KL for pair of Poisson(alpha) tilted by e^{+-u/2}
klp=lambda l1,l0: l1*sp.log(l1/l0)-l1+l0
KL=sp.simplify(klp(a*sp.exp(u/2),a)+klp(a*sp.exp(-u/2),a))
out['KL_pair']=str(sp.simplify(KL - (u*a*sp.sinh(u/2)-2*a*(sp.cosh(u/2)-1))).rewrite(sp.exp).simplify())
out['mean_pair']=str(sp.simplify(((a*sp.exp(u/2)-a*sp.exp(-u/2))/2 - a*sp.sinh(u/2)).rewrite(sp.exp)))
# 4. Legendre transform V*(B) = 2B asinh(B/4) - 2 sqrt(B^2+16) + 8
B=sp.symbols('B',positive=True); tau=2*sp.asinh(B/4)
Vs=sp.simplify(tau*B-8*(sp.cosh(tau/2)-1)); out['Vstar']=str(sp.simplify(Vs-(2*B*sp.asinh(B/4)-2*sp.sqrt(B**2+16)+8)))
# 5. prime-power sieve to X
X=10**7
sv=np.ones(X+1,bool); sv[:2]=False
for p in range(2,int(X**.5)+1):
    if sv[p]: sv[p*p::p]=False
primes=np.nonzero(sv)[0]
pp=[]
for p in primes:
    q=int(p); lp=math.log(p)
    while q<=X: pp.append((q,lp)); q*=int(p)
pp.sort(); N=np.array([q for q,_ in pp],float); Lam=np.array([l for _,l in pp])
alpha=2*Lam/(N-1); A=np.cumsum(alpha)
out['A_through_5']=float(A[np.searchsorted(N,5)]); out['A_through_7']=float(A[np.searchsorted(N,7)])
out['A_minus_2logx']=[(int(x),float(A[np.searchsorted(N,x,side='right')-1]-2*math.log(x))) for x in [1e2,1e3,1e4,1e5,1e6,1e7]]
# 6. Grok hypothesis test: on interval I_k=[log n_k, log n_{k+1}), (B,C) constant; inf_t (V(t)-tB) = -V*(B) at tau*=2asinh(B/4).
w=Lam/np.sqrt(N); Bk=np.cumsum(w); Ck=np.cumsum(w*np.log(N))
Vstar=lambda b: 2*b*np.arcsinh(b/4)-2*np.sqrt(b*b+16)+8
gap=Ck-Vstar(Bk)               # C_k - V*(B_k): Legendre-optimal coupling feasible on I_k iff >=0
tau=2*np.arcsinh(Bk/4); lo=np.log(N); hi=np.append(np.log(N[1:]),math.log(X+1))
inside=(tau>=lo)&(tau<hi)
V=lambda tt: 8*(np.cosh(tt/2)-1)
Tlo=V(lo)-lo*Bk+Ck; Thi=V(hi)-hi*Bk+Ck
Tmin_int=np.where(inside,gap,np.minimum(Tlo,Thi))   # exact min of convex T on I_k
out['n_intervals']=int(len(N)); out['frac_tau_inside']=float(inside.mean())
i=int(np.argmin(gap)); out['min_gap']={'value':float(gap[i]),'n_k':int(N[i]),'B':float(Bk[i]),'C':float(Ck[i])}
out['num_gap_negative']=int((gap<0).sum())
j=int(np.argmin(Tmin_int)); out['min_T']={'value':float(Tmin_int[j]),'n_k':int(N[j])}
out['num_T_negative']=int((Tmin_int<0).sum())
# relative size: gap / e^{t/2} on decades
out['gap_scaled']=[(int(x),float(gap[np.searchsorted(N,x,side='right')-1]/math.sqrt(x))) for x in [1e2,1e3,1e4,1e5,1e6,1e7]]
neg=np.nonzero(gap<0)[0][:15]; out['first_negative_gaps']=[(int(N[k]),float(gap[k]),bool(inside[k])) for k in neg]
zp=float(mp.diff(mp.zeta,0.5)/mp.zeta(0.5)); out['zeta_log_deriv_half']=zp
out['gap_vs_drift']=[(int(x),float(gap[np.searchsorted(N,x,side='right')-1]),zp*math.log(x),float(Tmin_int[np.searchsorted(N,x,side='right')-1])) for x in [1e2,1e3,1e4,1e5,1e6,1e7]]
dd=gap-zp*np.log(N); out['gap_minus_drift_range']=[float(dd[100:].min()),float(dd[100:].max())]
json.dump(out,open('h42b_poisson.json','w'),indent=1,default=str); print(json.dumps(out,indent=1,default=str))
