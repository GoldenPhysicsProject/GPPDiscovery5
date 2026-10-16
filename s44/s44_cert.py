"""S44 outer-bound certificates (Arb): for each sector take the Ritz extreme generalized eigenvector x of T2 v = mu W_zeta v,
round to rationals, and compute in Arb  r(x) = W_zeta(x)/T2(x)  with all entries from the exact S14 assembly.
 If T2(x)>0 then W(c)(x) = W_zeta(x) - (c-c_true)T2(x) < 0 for all c > c_true + r(x): CERTIFIED upper bound on the true upper gap.
 If T2(x)<0 then W(c)(x) < 0 for all c < c_true + r(x) (r<0): CERTIFIED upper bound on the true lower gap |r|.
These are outer bounds only (true window subset of [c_true - |r_-|, c_true + r_+]); no inner bound is claimed here."""
import sys, json; sys.path.insert(0,'.')
from s44_scan import *
def cert(M,h,prec):
    A,T,ns=build_all(M,h,prec)
    W=arb_mat(M,M)
    for i in range(M):
        for j in range(M):
            W[i,j]=A[i,j]-sum((2*arb(vonmangoldt(n)).log()/arb(n).sqrt())*T[n][i,j] for n in ns)
    row={"M":M,"h":h,"prime_powers":ns}
    for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
        w=tomp(W,I); t=tomp(T[2],I)
        Lc=mp.cholesky(w); Li=mp.inverse(Lc); S=Li*t*Li.T; S=(S+S.T)/2
        mu,V=mp.eigsy(S); out={}
        for side,idx in (("up",max(range(len(I)),key=lambda i:mu[i])),("down",min(range(len(I)),key=lambda i:mu[i]))):
            y=Li.T*V[:,idx]; sc=max(abs(v) for v in y); y=[v/sc for v in y]
            xs=[fmpq(int(mp.nint(v*mp.mpf(10)**70)),10**70) for v in y]
            xa=[arb(q.p)/q.q for q in xs]
            wv=sum(xa[a]*xa[b]*W[I[a],I[b]] for a in range(len(I)) for b in range(len(I)))
            tv=sum(xa[a]*xa[b]*T[2][I[a],I[b]] for a in range(len(I)) for b in range(len(I)))
            r=wv/tv
            out[side]={"T2x_sign":"+" if tv>0 else ("-" if tv<0 else "?"),"r_cert":r.str(12),"ritz_gap":mp.nstr(1/mu[idx],12)}
        row[nm]=out
    return row
if __name__=="__main__":
    M=int(sys.argv[1])
    for h in sys.argv[2].split(','): print(json.dumps(cert(M,h,int(sys.argv[4]) if len(sys.argv)>4 else 768)),flush=True)
