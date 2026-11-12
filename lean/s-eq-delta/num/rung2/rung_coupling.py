"""NUMERICAL (float64): RungStepCoupling written out at the rung 1 -> 2 step (n = 3 enters) and the
rung 2 -> 3 step (n = 4 enters).  No zeros used.

Old window V_k = L^2[-h_k, h_k] (in-tree rung window L = 2 h_k); new window [-h_{k+1}, h_{k+1}].
Basis: orthonormal Legendre modes on [-h_k, h_k] (V_k) and on each annulus piece [h_k, h_{k+1}],
[-h_{k+1}, -h_k] (exactly V_k^perp).  Block form Q = [[A, B], [B^T, D]].
(a) Galerkin of Zhu's lower form R_T (R_T <= Q, R_T increasing to Q in T) for T in (800, 3200).
(b) TRUE-form bilinear values Q(phi_1, psi), Q(psi, psi), split into pole / archimedean / prime
    pieces, in the time domain with the in-tree definitions (cross-correlation, Levy split).
Reported: lambda_k = lambda_min(A) and its mode phi_k; lambda_{k+1} = lambda_min(Q);
lambda_min of the Schur complement S = D - B^T A^{-1} B; for psi* = ground state of S and for the
worst psi for the 2x2 test (max_psi b^2/(lambda_k d) = (B^T phi)^T D^{-1} (B^T phi)/lambda_k):
b = Q(phi_k, psi), d = Q(psi, psi), ratio b^2/(lambda_k d)."""
import numpy as np, math, json, sys
from scipy.special import spherical_jn, psi as digamma, eval_legendre
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss
gamma=0.5772156649015329; Cinf=-gamma-math.log(math.pi)-math.pi/2-3*math.log(2)
def mangoldt(n):
    for p in range(2,n+1):
        if n%p==0:
            m=n
            while m%p==0: m//=p
            return (math.log(p) if m==1 else 0.0)
def comb(h):  # prime powers with log n < 2h
    return [(n,mangoldt(n)) for n in range(2,int(math.exp(2*h))+2) if math.log(n)<2*h-1e-15 and mangoldt(n)>0]
def basis(hk,hk1,K0,K1):
    # list of (c, w, k): piece centre, half-width, degree
    B=[(0.0,hk,k) for k in range(K0)]
    cA=(hk+hk1)/2; wA=(hk1-hk)/2
    B+= [(cA,wA,k) for k in range(K1)] + [(-cA,wA,k) for k in range(K1)]
    return B
def ev_mode(m,x):
    c,w,k=m; z=(x-c)/w; r=np.zeros_like(x,dtype=float); s=np.abs(z)<=1
    r[s]=math.sqrt((2*k+1)/(2*w))*eval_legendre(k,z[s]); return r
def RT_matrix(B,h1,T,terms):
    A=sum(2*l/math.sqrt(n) for n,l in terms); beta=math.log(T/(2*math.pi))-1/T-A
    g20,w20=leggauss(20); npan=int(T/0.5); a=np.arange(npan)*0.5
    t=(a[:,None]+(g20[None,:]+1)*0.25).ravel(); tw=np.tile(w20*0.25,npan)
    Psi=digamma(0.25+0.5j*t).real-math.log(math.pi)-sum(2*l/math.sqrt(n)*np.cos(t*math.log(n)) for n,l in terms)
    F=np.array([np.exp(-1j*t*c)*math.sqrt((2*k+1)/(2*w))*w*2*((-1j)**k)*spherical_jn(k,w*t) for c,w,k in B])
    M=np.real((np.conj(F)*(tw*(Psi-beta))) @ F.T)/math.pi
    xg,wg=leggauss(60)
    def proj(fun):
        out=[]
        for c,w,k in B:
            x=c+w*xg; out.append(float(np.sum(w*wg*ev_mode((c,w,k),x)*fun(x))))
        return np.array(out)
    ch=proj(lambda x: np.cosh(x/2)); sh=proj(lambda x: np.sinh(x/2))
    R=2*np.outer(ch,ch)-2*np.outer(sh,sh)+M+beta*np.eye(len(B))
    return (R+R.T)/2
# ---- true-form bilinear in the time domain ----
XG,WG=leggauss(40)
def make_f(B,coef):
    def f(x):
        x=np.atleast_1d(np.asarray(x,dtype=float)); r=np.zeros_like(x)
        for cf,m in zip(coef,B):
            if cf!=0.0: r+=cf*ev_mode(m,x)
        return r
    return f
def xcorr(f1,f2,u,bps):
    # int f1(y) f2(y-u) dy, piecewise Gauss on breakpoints bps and bps+u
    pts=sorted(set([p for p in bps]+[p+u for p in bps]))
    s=0.0
    for a,b in zip(pts[:-1],pts[1:]):
        if b-a<1e-15: continue
        y=(XG+1)*(b-a)/2+a; s+=float(np.sum(WG*(b-a)/2*f1(y)*f2(y-u)))
    return s
def Qbil(f1,f2,H,bps,terms):
    gs=lambda u: 0.5*(xcorr(f1,f2,u,bps)+xcorr(f1,f2,-u,bps))
    g0=gs(0.0)
    xg,wg=leggauss(200)
    def mom(f,sg):
        s=0.0
        for a,b in zip(bps[:-1],bps[1:]):
            y=(xg+1)*(b-a)/2+a; s+=float(np.sum(wg*(b-a)/2*f(y)*np.exp(sg*y/2)))
        return s
    pole=mom(f1,1)*mom(f2,-1)+mom(f1,-1)*mom(f2,1)
    kinks=sorted(set(abs(p-q) for p in bps for q in bps if 0<abs(p-q)<2*H))
    I1=quad(lambda u: math.exp(-u/2)*(g0-gs(u))/(1-math.exp(-2*u)),0,2*H,points=kinks,limit=800,epsabs=1e-16,epsrel=1e-12)[0]
    I2=g0*quad(lambda u: math.exp(-u/2)/(1-math.exp(-2*u)),2*H,80,epsabs=1e-17,epsrel=1e-14)[0]
    arch=Cinf*g0+2*(I1+I2)
    primes={n: l/math.sqrt(n)*2*gs(math.log(n)) for n,l in terms}
    return dict(total=pole+arch-sum(primes.values()),pole=pole,arch=arch,**{f"prime_{n}":v for n,v in primes.items()})
def step(name,hk,hk1,K0,K1,Ts=(800.0,3200.0)):
    B=basis(hk,hk1,K0,K1); n0=K0
    terms=comb(hk1); res=dict(name=name,h_old=hk,h_new=hk1,comb=[n for n,_ in terms],K0=K0,K1=K1)
    for T in Ts:
        R=RT_matrix(B,hk,T,terms)
        A=R[:n0,:n0]; Bm=R[:n0,n0:]; D=R[n0:,n0:]
        la,Va=np.linalg.eigh(A); lam_k=la[0]; phi=Va[:,0]
        lam_new=np.linalg.eigvalsh(R)[0]
        S=D-Bm.T@np.linalg.solve(A,Bm); ls,Vs=np.linalg.eigh(S)
        psi_s=Vs[:,0]; b_s=phi@Bm@psi_s; d_s=psi_s@D@psi_s
        bv=Bm.T@phi; y=np.linalg.solve(D,bv); psi_w=y/np.linalg.norm(y)
        b_w=phi@Bm@psi_w; d_w=psi_w@D@psi_w
        r=dict(T=T,lambda_old=float(lam_k),lambda_new=float(lam_new),schur_min=float(ls[0]),
               psi_schur=dict(b=float(b_s),d=float(d_s),ratio=float(b_s**2/(lam_k*d_s))),
               psi_worst=dict(b=float(b_w),d=float(d_w),ratio=float(b_w**2/(lam_k*d_w))),
               D_min=float(np.linalg.eigvalsh(D)[0]))
        print(name,"T=",T,json.dumps(r),flush=True); res[f"T{int(T)}"]=r
    # true-form values at T=3200 vectors
    phi_full=np.concatenate([phi,np.zeros(len(B)-n0)])
    bps=sorted(set([-hk1,-hk,hk,hk1]))
    fphi=make_f(B,phi_full)
    out={}
    for lab,ps in (("psi_schur",psi_s),("psi_worst",psi_w)):
        fpsi=make_f(B,np.concatenate([np.zeros(n0),ps]))
        bb=Qbil(fphi,fpsi,hk1,bps,terms); dd=Qbil(fpsi,fpsi,hk1,bps,terms)
        out[lab]=dict(b_parts=bb,d_parts=dd)
    aa=Qbil(fphi,fphi,hk1,bps,terms); out["Q_phi"]=aa
    for lab in ("psi_schur","psi_worst"):
        b=out[lab]["b_parts"]["total"]; d=out[lab]["d_parts"]["total"]
        out[lab]["ratio_true"]=b*b/(aa["total"]*d)
    print(name,"TRUE-form:",json.dumps(out),flush=True)
    res["true_form"]=out
    xs=np.linspace(-hk,hk,11); prof=fphi(xs); res["phi_profile"]=(prof/np.max(np.abs(prof))).round(4).tolist()
    return res
if __name__=="__main__":
    L2,L3,L4,L5=(math.log(n) for n in (2,3,4,5))
    out=[step("rung1->2 (n=3)",L3/2,L4/2,int(sys.argv[1]) if len(sys.argv)>1 else 30,int(sys.argv[3]) if len(sys.argv)>3 else 12),
         step("rung2->3 (n=4)",L4/2,L5/2,int(sys.argv[1]) if len(sys.argv)>1 else 30,int(sys.argv[3]) if len(sys.argv)>3 else 12)]
    json.dump(out,open(sys.argv[2] if len(sys.argv)>2 else "rung_coupling.json","w"),indent=1)
