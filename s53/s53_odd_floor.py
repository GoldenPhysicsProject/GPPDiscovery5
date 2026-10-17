"""Grok S53 H53b: odd-floor mechanism. S14/S44 exact Arb assembly (Legendre*(1-u^2)^4, M functions), split into
Pole, Arch, and each prime-power term -w_n T_n. For each sector, the minimizer v of W_zeta (Ritz, mpmath eigsy) and its
Rayleigh quotient decomposed: lam = v'Pv + v'Av - sum_n w_n v'T_n v  (v'v=1). NUMERICAL (assembly in Arb; eig in mpmath)."""
import sys, json, time; sys.path.insert(0,'.')
from s14_exact import *
import mpmath as mp
from flint import arb, fmpq
mp.mp.dps=int(sys.argv[3]) if len(sys.argv)>3 else 130
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
    A=arb_mat(M,M); P=arb_mat(M,M); T={n:arb_mat(M,M) for n in ns}
    for (j,k),g in Gs.items():
        co=[arb(c.p)/c.q for c in g.coeffs()]
        G0=co[0]; q=[-c for c in co[1:]]
        pole=4*ha*ha*sum(co[m]*m3[m] for m in range(len(co)))
        arch=ha*c0*G0 + ha*sum(q[m]*m1[m] for m in range(len(q))) - 2*ha*ha*sum(co[m]*m2[m] for m in range(len(co)))
        A[j,k]=arch; P[j,k]=pole
        for n in ns:
            v=arb(n).log()/ha
            T[n][j,k]=ha*sum(co[m]*v**m for m in range(len(co)))
    return A,P,T,ns
def tomp(Ma,I): return mp.matrix([[mp.mpf(Ma[i,j].mid().str(mp.mp.dps+20,radius=False)) for j in I] for i in I])
def weight(n): p=vonmangoldt(n); return 2*mp.log(p)/mp.sqrt(n)
def psd(X):
    try: mp.cholesky(X); return True
    except Exception: return False

def run(M,h,prec):
    t0=time.time(); A,P,T,ns=build_all(M,h,prec)
    row={"M":M,"h":h,"prime_powers":ns}
    for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
        a=tomp(A,I); p=tomp(P,I); Tm={n:tomp(T[n],I) for n in ns}
        W=a+p
        for n in ns: W-=weight(n)*Tm[n]
        ev,Q=mp.eigsy(W); i0=min(range(len(I)),key=lambda i: ev[i]); v=Q[:,i0]
        q=lambda X: (v.T*X*v)[0]
        dec={"pole":mp.nstr(q(p),8),"arch":mp.nstr(q(a),8)}
        for n in ns: dec["n=%d"%n]=mp.nstr(-weight(n)*q(Tm[n]),8)
        dec["sum_primes"]=mp.nstr(sum(-weight(n)*q(Tm[n]) for n in ns),8)
        evs=sorted(ev)
        row[nm]={"lam_min":mp.nstr(evs[0],10),"lam_2":mp.nstr(evs[1],10),"decomp_of_minimizer":dec,
                 "coef_mass_first4":mp.nstr(sum(v[i]**2 for i in range(4)),6)}
    row["wall_s"]=round(time.time()-t0,1); return row
if __name__=="__main__":
    M=int(sys.argv[1]); hs=sys.argv[2].split(","); prec=int(sys.argv[4]) if len(sys.argv)>4 else 1024
    for h in hs:
        print(json.dumps(run(M,mp.mpf(eval(h)) if False else h,prec)),flush=True)
