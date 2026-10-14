"""S13: coupling of prime 2 to the pole+arch negative directions. NUMERICAL (float Ritz, M=14), exact window fact noted.
For h < (1/2)log 3 the only prime power with log n < 2h is n=2, so W|C_h = A(h) - c*T2(h) EXACTLY, where
A = pole+arch, T2 = Fmat(log 2) (compressed translation by log 2, symmetrized), c_zeta = 2 log2/sqrt2 = sqrt2 log2.
F_theta: A_theta (shifted pole+arch), c_theta = (2^th+2^-th) log2/sqrt2.
Computes: admissible window [c_lo,c_hi] with A - c T2 >= 0 on the Ritz space, alignment of A's negative vectors with T2."""
import numpy as np, sys
from scipy.linalg import eigh
src=open("experiments/weil_small_support/weil_ritz.py").read().split("cases = [")[0]
sys.argv=[sys.argv[0],"14","0.3","0.3","1"]; ns={"__file__":"experiments/weil_small_support/weil_ritz.py"}; exec(src,ns)
M=ns["M"]; L2=np.log(2)
def T2(h):
    xq,wq=ns["xq0"]*h,ns["wq0"]*h; B=ns["basis"](xq,h); A=(B*wq)@ns["basis"](xq-L2,h).T; return (A+A.T)/2
def window(A,T,G):
    # c-range with A - cT >= 0 : generalized eigen of pencil; scan
    cs=np.linspace(-6,6,24001); ok=[c for c in cs if eigh(A-c*T,G,eigvals_only=True)[0]>=0]
    return (min(ok),max(ok)) if ok else None
for th in (0.0,0.25):
    c=(2**th+2**-th)*L2/np.sqrt(2)
    for h in (0.40,0.45,0.50,0.54):
        Q,G=ns["form"](h,"Z",th); T=T2(h); A=Q+c*T
        for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
            a,t,g=A[np.ix_(I,I)],T[np.ix_(I,I)],G[np.ix_(I,I)]
            wa,va=eigh(a,g); wt,vt=eigh(t,g)
            x=va[:,0]; y=vt[:,0]; ov=abs(x@g@y)/np.sqrt((x@g@x)*(y@g@y))
            nneg=int((wa<0).sum()); W=window(a,t,g)
            print(f"th={th} h={h:.2f} {nm}: #neg(A)={nneg} lamA={wa[0]:+.3e} <x,T2x>={x@t@x:+.3e} T2min={wt[0]:+.3e} overlap(botA,botT2)={ov:.3f} c={c:.4f} window={None if W is None else (round(W[0],4),round(W[1],4))} lamW={eigh(a-c*t,g,eigvals_only=True)[0]:+.3e}",flush=True)
