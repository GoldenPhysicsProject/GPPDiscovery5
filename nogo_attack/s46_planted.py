"""S46 H46a: planted off-line zero detection-support map on the S14/S44 exact assembly.
Weil form (S14 exact Arb, Legendre*(1-u^2)^4 basis, g supported in [-L,L], F=g*g~ in [-2L,2L]):
  W_zeta(F) = Pole + Arch - sum_n 2Lambda(n)/sqrt(n) F(log n) = sum_rho Fhat(gamma_rho),  Fhat(z)=int F(x)e^{izx}dx.
Planted model: replace on-line pair gamma_a, gamma_b (consecutive zeta zeros, both +-) by the quadruple 1/2+-delta+-i*gamma:
  Delta(F) = 4 int F cos(gamma x)cosh(delta x)dx - 2 int F [cos(gamma_a x)+cos(gamma_b x)] dx,
  Q_fake = W_zeta + Delta.  Entries: int F_jk phi = 2L^2 int_0^2 G_jk(v) phi(Lv) dv (F_jk even), moments
  I_m(a)=int_0^2 v^m e^{av}dv by the exact recursion I_m=(2^m e^{2a}-m I_{m-1})/a in Arb (rigorous balls).
Zeros gamma_a,gamma_b from GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt (>=30 digits), entered as balls of radius 1e-28."""
import sys, json, time; sys.path.insert(0,'.')
from s44_scan import *
from flint import acb, arb, arb_mat, fmpq, ctx
ZF='/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt'
HPZ={}
def zeros(N, hp=(1,2,10,11,23,24,50,51)):
    """table zeros; indices in hp replaced by mpmath zetazero at 70 digits (|err|<1e-60, cross-checked vs table to 1e-31)."""
    out=[]
    for i,l in enumerate(open(ZF)):
        if i>=N: break
        out.append(l.split()[1])
    import mpmath as _m
    for n in hp:
        if n<=N:
            if n not in HPZ:
                with _m.workdps(75): HPZ[n]=_m.nstr(_m.zetazero(n).imag,70)
            out[n-1]=HPZ[n]
    return out
def Gpolys(M):
    Gs={}
    for j in range(M):
        for k in range(j,M):
            if (j+k)%2==0: g=(G(j,k)+G(k,j))/2; Gs[j,k]=Gs[k,j]=g
    return Gs
def moms(a,D):
    e=(2*a).exp(); I=[(e-1)/a]; p2=acb(1)
    for m in range(1,D+1):
        p2*=2; I.append((p2*e-m*I[-1])/a)
    return I
ZEPS='1e-55'
def zball(s): return arb(s)   # exact table midpoint; table error |eps|<=1e-28 added via add_zero_radius

def delta_mat(M,h,Gs,terms,prec=1536,zerr=ZEPS):
    """terms: list of (coef arb, alpha acb) with phi(x)=Re sum coef*e^{alpha x}. Returns arb_mat of int F_jk phi."""
    old=ctx.prec; ctx.prec=prec; ha=arb(h)
    D=max(g.degree() for g in Gs.values())+1
    mu=[arb(0)]*(D+1)
    for c,al in terms:
        I=moms(acb(ha)*al,D)
        mu=[mu[m]+(acb(c)*I[m]).real for m in range(D+1)]
    X=arb_mat(M,M)
    for (j,k),g in Gs.items():
        co=[arb(c.p)/c.q for c in g.coeffs()]
        X[j,k]=2*ha*ha*sum(co[m]*mu[m] for m in range(len(co)))
    if zerr is not None:   # rigorous inflation for zero-table error: |d/dgamma int F_jk c cos(gamma x)| <= |c| int|x||F_jk| <= |c| 2L (2L)^2
        r=arb(0); 
        for c,al in terms:
            if c<0: r+=abs(c)*16*ha**3*arb(zerr)
        rr=arb(0,r.upper())
        for (j,k) in Gs: X[j,k]=X[j,k]+rr
    ctx.prec=old
    return X
def planted_terms(gam,dl,ga,gb):
    I=acb(0,1)
    return [(arb(2),I*gam+dl),(arb(2),I*gam-dl),(arb(-2),I*ga),(arb(-2),I*gb)]
def Wzeta(M,h,prec):
    A,T,ns=build_all(M,h,prec)
    W=arb_mat(M,M)
    for i in range(M):
        for j in range(M):
            W[i,j]=A[i,j]-sum((2*arb(vonmangoldt(n)).log()/arb(n).sqrt())*T[n][i,j] for n in ns)
    return W,ns
def secmin(X,I):
    x=tomp(X,I); ev,V=mp.eigsy(x); i=min(range(len(I)),key=lambda t:ev[t]); return ev[i],V[:,i]
def certify(X,I,v):
    sc=max(abs(t) for t in v); y=[t/sc for t in v]
    xs=[fmpq(int(mp.nint(t*mp.mpf(10)**60)),10**60) for t in y]; xa=[arb(q.p)/q.q for q in xs]
    val=sum(xa[a]*xa[b]*X[I[a],I[b]] for a in range(len(I)) for b in range(len(I)))
    return val
