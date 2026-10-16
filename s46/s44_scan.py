"""S44 H44a: prime-2 (k=1) coupling window vs support L (=h, test fns g supported in [-h,h], F=g*g~ in [-2h,2h]).
Weil form (S14 exact Arb assembly, Legendre*(1-u^2)^4 basis, M functions):
  W(c) = Pole + Arch - sum_{n>=2, log n<2h} w_n T_n  with w_n = 2 Lambda(n)/sqrt(n), T_n = F(log n),
  and the n=2 weight w_2 = sqrt2*log2 replaced by c.  W(c) = W_zeta - (c - c_true) T_2.
Per sector (even/odd), Ritz: generalized eigenproblem T_2 v = mu W_zeta v (W_zeta PD on Ritz space);
  delta_+ = 1/max(mu>0)  (upper gap c_+ - c_true), delta_- = 1/|min(mu<0)|  (lower gap).
Ritz space is a subspace, so Ritz window is an OUTER bound for the true window (gaps can only shrink as M grows).
Optional tower mode: 2^k weights scaled a^k (Euler factor (1-2^-s)/(1-a 2^-s)), bisection on a.  NUMERICAL."""
import sys, json, time; sys.path.insert(0,'.')
from s14_exact import *
import mpmath as mp
from flint import arb, fmpq
mp.mp.dps=int(sys.argv[3]) if len(sys.argv)>3 else 110
def vonmangoldt(n):
    for p in range(2,n+1):
        if n%p==0:
            m=n
            while m%p==0: m//=p
            return (p if m==1 else None)
    return None
def build_all(M,h,prec):
    ctx.prec=prec; ha=arb(h)
    Gs={}
    for j in range(M):
        for k in range(j,M):
            if (j+k)%2==0: g=(G(j,k)+G(k,j))/2; Gs[j,k]=Gs[k,j]=g
    D=max(g.degree() for g in Gs.values())+1
    def mom(phi):
        return [acb.integral(lambda z,an,m=m: z**m*phi(z),0,2).real for m in range(D+1)]
    hc=acb(ha)
    m1=mom(lambda z: (-hc*z).exp()/ (acb(0,1)*hc*z).sinc() )
    m2=mom(lambda z: (-hc*z/4).exp()*(acb(0,1)*hc*z*3/4).sinc()*fmpq(3,4)/(acb(0,1)*hc*z).sinc())
    m3=mom(lambda z: (hc*z/2).cosh())
    c0=-arb.const_euler()-arb.pi().log()-(1-(-4*ha).exp()).log()
    ns=[n for n in range(2,int(mp.exp(2*mp.mpf(h)))+2) if vonmangoldt(n) and mp.log(n)<2*mp.mpf(h)]
    A=arb_mat(M,M); T={n:arb_mat(M,M) for n in ns}
    for (j,k),g in Gs.items():
        co=[arb(c.p)/c.q for c in g.coeffs()]
        G0=co[0]; q=[-c for c in co[1:]]
        pole=4*ha*ha*sum(co[m]*m3[m] for m in range(len(co)))
        arch=ha*c0*G0 + ha*sum(q[m]*m1[m] for m in range(len(q))) - 2*ha*ha*sum(co[m]*m2[m] for m in range(len(co)))
        A[j,k]=pole+arch
        for n in ns:
            v=arb(n).log()/ha
            T[n][j,k]=ha*sum(co[m]*v**m for m in range(len(co)))
    return A,T,ns
def tomp(Ma,I): return mp.matrix([[mp.mpf(Ma[i,j].mid().str(mp.mp.dps+20,radius=False)) for j in I] for i in I])
def weight(n): p=vonmangoldt(n); return 2*mp.log(p)/mp.sqrt(n)
def psd(X):
    try: mp.cholesky(X); return True
    except Exception: return False
def run(M,h,prec,tower=False):
    t0=time.time(); A,T,ns=build_all(M,h,prec)
    maxrad=max(float(A[i,j].rad()) for i in range(M) for j in range(M))
    row={"M":M,"h":h,"prime_powers":ns,"maxrad_A":maxrad,"c_true":mp.nstr(weight(2),20)}
    for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
        a=tomp(A,I); Tm={n:tomp(T[n],I) for n in ns}
        W=a.copy()
        for n in ns: W-=weight(n)*Tm[n]
        ev=mp.eigsy(W,eigvals_only=True); lmin=min(ev)
        r={"lam_min_Wzeta":mp.nstr(lmin,8)}
        if lmin<=0: r["note"]="W_zeta NOT PD on Ritz space"; row[nm]=r; continue
        Lc=mp.cholesky(W); Li=mp.inverse(Lc)
        S=Li*Tm[2]*Li.T; S=(S+S.T)/2
        mu=mp.eigsy(S,eigvals_only=True); mx,mn=max(mu),min(mu)
        r["mu_max"]=mp.nstr(mx,10); r["mu_min"]=mp.nstr(mn,10)
        r["gap_up"]=mp.nstr(1/mx,10) if mx>0 else "inf"
        r["gap_down"]=mp.nstr(-1/mn,10) if mn<0 else "inf"
        if tower:
            k2=[n for n in ns if vonmangoldt(n)==2]
            base=a.copy()
            for n in ns:
                if n not in k2: base-=weight(n)*Tm[n]
            def ok(aa): 
                X=base.copy()
                for n in k2:
                    k=int(round(mp.log(n)/mp.log(2))); X-=aa**k*(2*mp.log(2)/mp.sqrt(n))*Tm[n]
                return psd(X)
            ed={}
            for side,far in (("a_up",mp.mpf(1.5)),("a_down",mp.mpf(0.5))):
                lo,hi=mp.mpf(1),far
                if ok(hi): ed[side]="beyond "+mp.nstr(far,3); continue
                for _ in range(60):
                    m=(lo+hi)/2
                    if ok(m): lo=m
                    else: hi=m
                ed[side]=mp.nstr(lo-1,10)
            r["tower_2k"]={"k_in":[int(round(mp.log(n)/mp.log(2))) for n in k2],**ed}
        row[nm]=r
    row["sec"]=round(time.time()-t0,1)
    return row
if __name__=="__main__":
    M=int(sys.argv[1]); hs=sys.argv[2].split(','); prec=int(sys.argv[4]) if len(sys.argv)>4 else 768
    tower=len(sys.argv)>5 and sys.argv[5]=="tower"
    for h in hs: print(json.dumps(run(M,h,prec,tower)),flush=True)
