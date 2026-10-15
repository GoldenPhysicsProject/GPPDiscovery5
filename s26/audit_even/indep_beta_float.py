"""Grok S26 independent float audit of Astra even h=2/5 candidate.
Independent formulas: per-segment (non-telescoped) PL cosine transform, hand-derived
closed form for int_0^h cosh(t/2)cos(xi t), Gauss-Legendre panels, scipy complex digamma.
Not a certificate; cross-check of beta_r, objective and margin. RH OPEN."""
import json, sys, numpy as np
from fractions import Fraction as F
from scipy.special import digamma
from numpy.polynomial.legendre import leggauss
d=json.load(open(sys.argv[1])); e=d['extension']
h=d['h_q'][0]/d['h_q'][1]; H=e['H_q'][0]/e['H_q'][1]; N=e['N']
lam=d['lambda_q'][0]/d['lambda_q'][1]; den=d['den']
C=d['C_num']/den; a=np.array(d['a_num'],float)/den; K=d['K']
D=d['D_q'][0]/d['D_q'][1]; Dd=d['Dd_q'][0]/d['Dd_q'][1]; c=d['c'][0]/d['c'][1]
ys=np.array(e['alpha_num'],float)/e['alpha_den']; ys=np.append(ys,0.0)
zs=np.array(e['z_num'],float)/e['z_den']
assert len(zs)==len(ys)==N+2
u=2*h+D*np.arange(K); kap=np.log(2)/np.sqrt(2)
def m(x):
    x=np.asarray(x,float)
    w=digamma(0.25+0.5j*x).real-np.log(np.pi)-2*kap*np.cos(x*np.log(2))
    Aw=2*h+Dd/2; b=Dd/2
    with np.errstate(all='ignore'):
        s1=np.where(x==0,Aw,np.sin(Aw*x)/np.where(x==0,1,x)); s2=np.where(x==0,1,np.sin(b*x)/np.where(x==0,1,b*x))
    w+=C*2*s1*s2
    w-=np.cos(np.outer(x,u))@a
    return w
def S(x):
    x=np.asarray(x,float); out=np.zeros_like(x)
    nz=x!=0; xx=x[nz]
    # interior: int_0^h (1 - lam cosh(t/2)) cos(xt)
    I=np.sin(h*xx)/xx - lam*(0.5*np.sinh(h/2)*np.cos(h*xx)+xx*np.cosh(h/2)*np.sin(h*xx))/(0.25+xx**2)
    # exterior per segment: int (y0+s(t-z0))cos = [(y0+s(t-z0)) sin/x + s cos/x^2]
    for i in range(N+1):
        z0,z1=zs[i],zs[i+1]; y0,y1=ys[i],ys[i+1]; s=(y1-y0)/(z1-z0)
        I+= (y1*np.sin(xx*z1)-y0*np.sin(xx*z0))/xx + s*(np.cos(xx*z1)-np.cos(xx*z0))/xx**2
    out[nz]=I
    out[~nz]=h-2*lam*np.sinh(h/2)+np.sum((ys[:-1]+ys[1:])*(zs[1:]-zs[:-1])/2)
    return out
# spot check S vs brute Gauss quadrature
gx,gw=leggauss(40)
def S_brute(x):
    tot=0.0
    pts=np.concatenate([np.linspace(0,h,201),zs[1:]])
    for i in range(len(pts)-1):
        A_,B_=pts[i],pts[i+1]; n=max(1,int(np.ceil(abs(x)*(B_-A_)/3)))
        for k in range(n):
            a0=A_+(B_-A_)*k/n; b0=A_+(B_-A_)*(k+1)/n
            t=(a0+b0)/2+(b0-a0)/2*gx
            r=np.where(t<h,1-lam*np.cosh(t/2),np.interp(t,zs,ys))
            tot+=np.sum(gw*r*np.cos(x*t))*(b0-a0)/2
    return tot
spot={x:(float(S(np.array([x]))[0]),S_brute(x)) for x in [0.0,0.37,2.0,11.0,64.5,250.0,499.0]}
# integrate S^2/m on [0,XB] with GL panels
XB=float(sys.argv[2]) if len(sys.argv)>2 else 500.0
pw=0.02; gx2,gw2=leggauss(16)
edges=np.arange(0,XB+1e-12,pw); bint=0.0; s2int=0.0; mmin=np.inf
for j in range(0,len(edges)-1,2000):
    L=edges[j:min(j+2001,len(edges))]
    mids=(L[:-1]+L[1:])/2; hw=(L[1:]-L[:-1])/2
    xs=(mids[:,None]+hw[:,None]*gx2[None,:]).ravel(); ws=(hw[:,None]*gw2[None,:]).ravel()
    mv=m(xs); sv=S(xs); mmin=min(mmin,mv.min())
    bint+=np.sum(ws*sv**2/mv); s2int+=np.sum(ws*sv**2)
vnorm2=h-4*lam*np.sinh(h/2)+lam**2*(h+np.sinh(h))/2 + np.sum((zs[1:]-zs[:-1])*(ys[:-1]**2+ys[:-1]*ys[1:]+ys[1:]**2)/3)
tailS2=np.pi/2*vnorm2-s2int
# tail estimate: m ~ m(XB) at least; use asymptotic m~log(x/2)-log pi for estimate
mXB=m(np.array([XB]))[0]
beta=4/np.pi*(bint+tailS2/mXB); beta_lo=4/np.pi*bint
obj=beta+lam**2/2; margin=1-C*obj
res=dict(h=h,H=H,N=N,lam=lam,C=C,c=c,XB=XB,beta_float_upperish=beta,beta_without_tail=beta_lo,tail_S2=tailS2,
 vnorm2=vnorm2,objective=obj,margin=margin,coercivity=margin*c,min_m_grid=float(mmin),m_minus_c_min=float(mmin-c),
 spot_S_formula_vs_brute={str(k):v for k,v in spot.items()},max_spot_abs_err=max(abs(a_-b_) for a_,b_ in spot.values()))
print(json.dumps(res,indent=1))
json.dump(res,open(f'indep_float_h04_XB{int(XB)}.json','w'),indent=1)
