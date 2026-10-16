"""Grok S42 H42a: independent audit of Codex spin-1/spin-2 mirrored CFS claims (PR11 6215940827, 6215318887).
sympy exact algebra + mpmath (dps 60) numerics. RH OPEN."""
import sympy as sp, mpmath as mp, json
out = {}
# --- 1. CFS Lagrangian normalization L = sum|l|^2 - (1/(2n)) (sum|l|)^2 = (1/(4n)) sum_{i,j}(|l_i|-|l_j|)^2
x1,x2 = sp.symbols('x1 x2', nonnegative=True)
n=1; L1 = x1**2+x2**2 - sp.Rational(1,2*n)*(x1+x2)**2
out['L_n1'] = str(sp.factor(sp.expand(L1)))           # (x1-x2)^2/2
n=2; ev=[x1,x1,x2,x2]; L2 = sum(e**2 for e in ev) - sp.Rational(1,2*n)*sum(ev)**2
out['L_n2_c1c1c2c2'] = str(sp.factor(sp.expand(L2)))  # (x1-x2)^2
ys=sp.symbols('y0:4',nonnegative=True)
out['pairwise_form_equal_n2'] = sp.simplify(sum(e**2 for e in ys)-sp.Rational(1,4)*sum(ys)**2 - sp.Rational(1,8)*sum((a-b)**2 for a in ys for b in ys))==0
# --- 2. J O J O^dagger = (|a|^2-|b|^2) I
ar,ai,br,bi=sp.symbols('ar ai br bi',real=True); a=ar+sp.I*ai; b=br+sp.I*bi
O=sp.Matrix([[a,b],[sp.conjugate(b),sp.conjugate(a)]]); J=sp.diag(1,-1)
out['JOJOdag_minus_scalar_zero'] = sp.simplify(sp.expand(J*O*J*O.H - (a*sp.conjugate(a)-b*sp.conjugate(b))*sp.eye(2)))==sp.zeros(2)
# also check O really is V_v^dag V_w for random vectors with coordinate conjugation (numeric, dim 6)
import random; random.seed(42)
v=[complex(random.gauss(0,1),random.gauss(0,1)) for _ in range(6)]; w=[complex(random.gauss(0,1),random.gauss(0,1)) for _ in range(6)]
ip=lambda p,q: sum(pi.conjugate()*qi for pi,qi in zip(p,q)); C=lambda p:[z.conjugate() for z in p]
Onum=[[ip(v,w),ip(v,C(w))],[ip(C(v),w),ip(C(v),C(w))]]
out['O_structure_numeric_err']=max(abs(Onum[1][0]-Onum[0][1].conjugate()),abs(Onum[1][1]-Onum[0][0].conjugate()))
# --- 3. Codex spin-1 bridge: (A_h - m I)^2 = h I, L=2max(h,0)
h=sp.symbols('h',real=True); m=1+sp.Abs(h); aa=(h+1)/2; vv=(h-1)/2
x=sp.diag(1,-1); y=sp.Matrix([[m+aa,vv],[vv,aa-m]]); A=x*y
out['spin1_bridge_(A-mI)^2-hI'] = [sp.simplify((A-m*sp.eye(2))**2 - h*sp.eye(2)).subs(h,hv)==sp.zeros(2) for hv in [-3,-sp.Rational(1,2),0,sp.Rational(1,3),5]]
out['spin1_bridge_detY'] = str(sp.simplify(y.det()))
Ls=[]
for hv in [-2.0,-0.1,0.0,0.3,4.0]:
    lam=sp.Matrix(A.subs(h,hv)).eigenvals(); ls=[complex(k) for k,mu in lam.items() for _ in range(mu)]
    Ls.append((hv, abs(ls[0])**2+abs(ls[1])**2-0.5*(abs(ls[0])+abs(ls[1]))**2, 2*max(hv,0)))
out['spin1_bridge_L_vs_2max(h,0)']=Ls
# --- 4. curvature identity Var_full - Var_sq = 4 (log zeta)''(2 sigma) = 4 sum_{p,k} k (log p)^2 p^(-2 k sigma)
mp.mp.dps=60
sig=mp.mpf(2)
lz=lambda s: mp.log(mp.zeta(s))
def d2logzeta(s):  # exact via zeta derivatives
    z=mp.zeta(s); z1=mp.zeta(s,derivative=1); z2=mp.zeta(s,derivative=2); return (z2*z-z1**2)/z**2
D=d2logzeta(2*sig)
Var1=d2logzeta(sig); Var2=mp.diff(lambda s: mp.log(mp.zeta(s)/mp.zeta(2*s)),sig,2)
out['Var1-Var2']=mp.nstr(Var1-Var2,40); out['4D']=mp.nstr(4*D,40); out['identity_err']=mp.nstr(abs(Var1-Var2-4*D),5)
# prime-sum form: k-weighted == sum Lambda(n) log n n^(-2sigma)
from sympy import primerange
N=200000; ps=list(primerange(2,N)); s2=mp.mpf(4); S=mp.mpf(0)
for p in ps:
    lp=mp.log(p); k=1; pk=mp.mpf(p)**(-s2)
    while pk>mp.mpf(10)**-70:
        S+=k*lp**2*pk; k+=1; pk*=mp.mpf(p)**(-s2)
tail=mp.mpf(1)/(N**2)  # crude: sum_{n>N} log^2 n n^-4 < log^2N/(3N^3)*~; we just quote
out['prime_sum_k_weighted_N2e5']=mp.nstr(S,30); out['prime_sum_minus_D']=mp.nstr(S-D,5)
# Lambda-weighted (WRONG variant) for contrast: sum_{p,k} (log p)^2 p^{-2ks}  (missing k)
S_wrong=sum(mp.log(p)**2*mp.mpf(p)**(-s2)/(1-mp.mpf(p)**(-s2)) for p in ps[:2000])
out['contrast_without_k']=mp.nstr(S_wrong,20)
# --- 5. self-action t^4 coefficient
def K(i,sg,u):
    if i==1: return mp.zeta(sg+1j*u)/mp.zeta(sg)
    return (mp.zeta(sg+1j*u)/mp.zeta(2*sg+2j*u))/(mp.zeta(sg)/mp.zeta(2*sg))
def L2(sg,t,s):
    c1=abs(K(1,sg,t-s))**2-abs(K(1,sg,t+s))**2; c2=abs(K(2,sg,t-s))**2-abs(K(2,sg,t+s))**2
    return (abs(c1)-abs(c2))**2, c1, c2
rows=[]
for tv in ['0.1','0.05','0.025','0.01','0.001']:
    t=mp.mpf(tv); Lv,c1,c2=L2(sig,t,t); pred=256*D**2*t**4
    rows.append({'t':tv,'L':mp.nstr(Lv,20),'pred':mp.nstr(pred,20),'ratio':mp.nstr(Lv/pred,12),'(1-ratio)/t^2':mp.nstr((1-Lv/pred)/t**2,8),
                 'c1/(4Var1 t^2)':mp.nstr(c1/(4*Var1*t**2),12),'c2/(4Var2 t^2)':mp.nstr(c2/(4*Var2*t**2),12)})
out['self_action']=rows
# --- 6. delta^-4 renormalized off-diagonal limit sigma=1+delta
def Alim(t,s):
    u=t-s; v=t+s; z2=mp.zeta(2)
    A1=abs(mp.zeta(1+1j*u))**2-abs(mp.zeta(1+1j*v))**2
    A2=z2**2*(abs(mp.zeta(1+1j*u)/mp.zeta(2+2j*u))**2-abs(mp.zeta(1+1j*v)/mp.zeta(2+2j*v))**2)
    return (abs(A1)-abs(A2))**2, A1, A2
lim=[]
for (t,s) in [(1,0.5),(2,0.7),(3,1)]:
    t=mp.mpf(t); s=mp.mpf(s); Lim,A1,A2=Alim(t,s); r={'t':float(t),'s':float(s),'limit':mp.nstr(Lim,15),'A1':mp.nstr(A1,10),'A2':mp.nstr(A2,10),'by_delta':[]}
    for dl in ['1e-2','1e-3','1e-4','1e-5']:
        d=mp.mpf(dl); Lv,_,_=L2(1+d,t,s); r['by_delta'].append({'delta':dl,'d^-4 L':mp.nstr(Lv/d**4,15),'relerr':mp.nstr(Lv/d**4/Lim-1,6),'relerr/delta':mp.nstr((Lv/d**4/Lim-1)/d,6)})
    lim.append(r)
out['delta_limit']=lim
json.dump(out,open('h42a_cfs_audit.json','w'),indent=1,default=str); print(json.dumps(out,indent=1,default=str))
