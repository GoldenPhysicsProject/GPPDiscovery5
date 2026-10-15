"""S24 [design, float, UNTRUSTED]: PL odd extension of v=sinh(x/2) beyond the window with a
general knot vector (graded at h, long range up to H ~ U+h), to minimise
beta = (1/2pi) int |v_w^|^2 / m. For odd functions v_w^ = -2i S(xi), S = int_0^inf v_w sin(xi x).
beta = (4/pi) int_0^inf S^2/m dxi (as in S23)."""
import numpy as np, scipy.special as sp, json, sys, time
l2=np.log(2); kap=l2/np.sqrt(2)
def S0(xi,h):
    a=0.5; num=a*np.cosh(a*h)*np.sin(xi*h)-xi*np.sinh(a*h)*np.cos(xi*h); return num/(a*a+xi*xi)
def m_fn(xi,d):
    h=float(d["h"]); a=np.array(d["a_num"],float)/d["den"]; C=d["C_num"]/d["den"]
    u=2*h+d["D"]*np.arange(d["K"]); Dd=d["Dd"]; A=2*h+Dd/2
    out=sp.psi(0.25+0.5j*xi).real-np.log(np.pi)-2*kap*np.cos(xi*l2)+C*2*A*np.sinc(A*xi/np.pi)*np.sinc(Dd*xi/(2*np.pi))
    for k in range(0,len(a),64): out-=np.cos(np.outer(xi,u[k:k+64]))@a[k:k+64]
    return out
def lin_sin(a,b,lo,hi,xi):
    f=lambda x: -(a*x+b)*np.cos(xi*x)/xi+a*np.sin(xi*x)/(xi*xi)
    return f(hi)-f(lo)
def tent(xi,l,mid,r):
    out=np.zeros_like(xi)
    if mid>l: out+=lin_sin(1/(mid-l),-l/(mid-l),l,mid,xi)
    if r>mid: out+=lin_sin(-1/(r-mid),r/(r-mid),mid,r,xi)
    return out
def knots(h,H,N,grade):
    t=np.linspace(0,1,N+2); return h+(H-h)*t**grade
def design(d,zs,dxi=0.01,Xi=1500.0,ridge=1e-12):
    h=float(d["h"]); xs=np.arange(dxi/2,Xi,dxi); mv=m_fn(xs,d); W=4/(np.pi*mv)*dxi
    base=S0(xs,h)+np.sinh(h/2)*tent(xs,zs[0],zs[0],zs[1])
    Phi=np.column_stack([tent(xs,zs[i-1],zs[i],zs[i+1]) for i in range(1,len(zs)-1)])
    A=(Phi*W[:,None]).T@Phi; b=(Phi*W[:,None]).T@base; c0=float(base*W@base)
    al=-np.linalg.solve(A+ridge*np.eye(len(b)),b)
    bc=c0+2*b@al+al@A@al
    # crude tail: S ~ O(1/xi^2) beyond Xi; m >= log
    return dict(beta=float(bc),two_beta=float(2*bc),alpha=al.tolist(),zs=list(map(float,zs)))
if __name__=="__main__":
    d=json.load(open(sys.argv[1])); h=float(d["h"]); res=[]
    for H in map(float,sys.argv[2].split(",")):
      for N in map(int,sys.argv[3].split(",")):
        for g in map(float,sys.argv[4].split(",")):
          t0=time.time(); o=design(d,knots(h,H,N,g)); o.update(H=H,N=N,grade=g)
          print(f"H={H} N={N} grade={g}: 2beta={o['two_beta']:.5f} ({time.time()-t0:.0f}s)",flush=True); res.append(o)
    best=min(res,key=lambda r:r['two_beta']); print("BEST",best['H'],best['N'],best['grade'],best['two_beta'])
    json.dump(best,open(sys.argv[5] if len(sys.argv)>5 else "best_design.json","w"))

def knots2(h,H,N1,N2,g1,g2):
    """graded at h (power g1), at h+log2 from both sides, then graded at h+log2 on the right (power g2)."""
    p=h+l2; t=np.linspace(0,1,N1+1)
    # two-sided grading on [h,p]: map t -> beta-like
    s=np.where(t<0.5,0.5*(2*t)**g1,1-0.5*(2*(1-t))**g1)
    z1=h+(p-h)*s
    t2=np.linspace(0,1,N2+1)[1:]; z2=p+(H-p)*t2**g2
    return np.r_[z1,z2]
def sweep2(d,H,N1s,N2s,g1s,g2s):
    res=[]
    for N1 in N1s:
     for N2 in N2s:
      for g1 in g1s:
       for g2 in g2s:
        zs=knots2(float(d["h"]),H,N1,N2,g1,g2); o=design(d,zs); o.update(H=H,N1=N1,N2=N2,g1=g1,g2=g2)
        print(f"H={H} N1={N1} N2={N2} g1={g1} g2={g2}: 2beta={o['two_beta']:.5f}",flush=True); res.append(o)
    return min(res,key=lambda r:r['two_beta'])

def design_tail(d,zs,dxi=0.01,Xi=1500.0,ridge=1e-12):
    """beta with honest tail: int_X^inf S^2/m <= (1/mX)[(pi/2)||v_w||^2 - int_0^X S^2]."""
    h=float(d["h"]); xs=np.arange(dxi/2,Xi,dxi); mv=m_fn(xs,d)
    xt=np.linspace(Xi,4*Xi,20001); mX=float(m_fn(xt,d).min())
    W=4/np.pi*(1/mv-1/mX)*dxi
    y0=np.sinh(h/2); n=len(zs)-2
    base=S0(xs,h)+y0*tent(xs,zs[0],zs[0],zs[1])
    Phi=np.column_stack([tent(xs,zs[i-1],zs[i],zs[i+1]) for i in range(1,n+1)])
    A=(Phi*W[:,None]).T@Phi; b=(Phi*W[:,None]).T@base; c0=float(base*W@base)
    # mass matrix of PL on knots with values (y0, y1..yn, 0)
    dz=np.diff(zs); Mfull=np.zeros((n+2,n+2))
    for i in range(n+1):
        Mfull[i,i]+=dz[i]/3; Mfull[i+1,i+1]+=dz[i]/3; Mfull[i,i+1]+=dz[i]/6; Mfull[i+1,i]+=dz[i]/6
    k=2/mX  # (4/pi)*(1/mX)*(pi/2)
    A+=k*Mfull[1:n+1,1:n+1]; b+=k*y0*Mfull[1:n+1,0]; c0+=k*(y0*y0*Mfull[0,0]+(np.sinh(h)-h)/2)
    al=-np.linalg.solve(A+ridge*np.eye(n),b); bc=c0+2*b@al+al@A@al
    return dict(beta=float(bc),two_beta=float(2*bc),alpha=al.tolist(),zs=list(map(float,zs)),mX=mX,Xi=Xi)
