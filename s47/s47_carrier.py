"""S47 H47b: carrier-modulated Ritz basis for the S46 planted-quartet model (NUMERICAL, float64 + zero-sum representation).
Q_fake(g) = sum_{zeta zeros n != a,b} 2|ghat(gamma_n)|^2 + tail + sum_{z in {+-gamma+-i delta}} ghat(z) ghat(-z),
  ghat(z) = int g(x) e^{izx} dx, g real supported in [-L,L].  (= W_zeta + Delta of S46 for g in the S14 span.)
Basis: p_j(u)=P_j(u)(1-u^2)^4, u=x/L: unmodulated j<K0, and p_j(u){cos,sin}(omega x) j<K, omega=gamma.
ghat_0j(z) = L sum_m c_jm 2 i^m j_m(zL) (Legendre expansion of p_j, spherical Bessel, complex z OK).
Tail beyond zero N: integral of (1/2pi)log(t/2pi) * 2|ghat|^2 (Riemann-von Mangoldt density)."""
import sys, json, numpy as np
from numpy.polynomial import legendre as Lg
from scipy.special import spherical_jn
from scipy.linalg import eigh
ZF='/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt'
NZ=int(__import__("os").environ.get("NZ","5000"))
GZ=np.array([float(l.split()[1]) for _,l in zip(range(NZ),open(ZF))])
W4=Lg.poly2leg([1,0,-2,0,1])**0  # placeholder
W4=Lg.poly2leg(np.polynomial.polynomial.polypow([1,0,-1],4))
def legco(j): return Lg.legmul(np.eye(j+1)[j],W4)
def ghat0(j,z,L):
    c=legco(j); z=np.asarray(z,dtype=complex); out=np.zeros(z.shape,complex)
    for m,cm in enumerate(c):
        if cm!=0: out+=cm*2*(1j**m)*spherical_jn(m,z*L)
    return L*out
def basis_hat(z,L,om,K,K0):
    rows=[ghat0(j,z,L) for j in range(K0)]
    for j in range(K):
        a,b=ghat0(j,z+om,L),ghat0(j,z-om,L)
        rows.append((a+b)/2); rows.append((a-b)/(2j))
    return np.array(rows)
def gram(L,om,K,K0,nq=4000):
    x,w=Lg.leggauss(nq); u=x; xx=L*x; ww=L*w
    B=[Lg.legval(u,legco(j)) for j in range(K0)]
    for j in range(K):
        p=Lg.legval(u,legco(j)); B+= [p*np.cos(om*xx), p*np.sin(om*xx)]
    B=np.array(B); return (B*ww)@B.T
def forms(L,om,K,K0,gam,dl,ia,ib):
    keep=np.ones(NZ,bool); keep[[ia-1,ib-1]]=False
    H=basis_hat(GZ[keep],L,om,K,K0)
    P=2*np.real(H@np.conj(H).T)
    # tail beyond gamma_20000~18047 omitted: positive, |ghat|^2 ~ t^-10 => contribution << 1e-20 (only could create false detections of that size; threshold is -1e-9)
    B=np.zeros_like(P)
    for z in (gam+1j*dl,gam-1j*dl,-gam+1j*dl,-gam-1j*dl):
        hp,hm=basis_hat(np.array([z]),L,om,K,K0)[:,0],basis_hat(np.array([-z]),L,om,K,K0)[:,0]
        X=np.outer(hp,hm); B+=np.real((X+X.T)/2)
    return P,B
def lmin(L,om,K,K0,gam,dl,ia,ib):
    P,B=forms(L,om,K,K0,gam,dl,ia,ib); G=gram(L,om,K,K0)
    s=1/np.sqrt(np.diag(G)); G2=G*np.outer(s,s); w,U=np.linalg.eigh(G2); k=w>1e-11*w.max()
    T=(U[:,k]/np.sqrt(w[k]))*s[:,None]   # orthonormal (L^2) coordinates on the numerically independent span
    ev0=np.linalg.eigvalsh(T.T@P@T); ev=np.linalg.eigvalsh(T.T@(P+B)@T)
    return float(ev[0]),float(ev0[0]),int(k.sum())
if __name__=="__main__":
    K=int(sys.argv[1]); K0=int(sys.argv[2])
    cfgs={"p1-2":(1,2,None),"p10-11":(10,11,None),"DHpt":(23,24,85.699),"p50-51":(50,51,None)}
    dls=[0.1,0.308]; Ls=[float(v) for v in sys.argv[3].split(',')]
    for nm,(ia,ib,g) in cfgs.items():
        gam=g if g else (GZ[ia-1]+GZ[ib-1])/2
        for dl in dls:
            row={"cfg":nm,"gamma":gam,"delta":dl,"K":K,"K0":K0,"omega":gam,"scan":[]}
            onset=None
            for L in Ls:
                l,l0,cn=lmin(L,gam,K,K0,gam,dl,ia,ib); row["scan"].append([L,l,l0])
                if l<-1e-9 and onset is None: onset=L
            row["onset_upper"]=onset; print(json.dumps(row),flush=True)
