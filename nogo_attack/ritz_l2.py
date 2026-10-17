"""Rigorous L2-normalised Ritz UPPER bounds for the zeta Weil form, per parity sector, on [-L,L].
Basis g_j(u)=p_j(u/L), p_j=P_j*(1-x^2)^4 (S14). Gram_jk = ||g||^2 pairing = L*G_jk(0) (F_jk(0)=int g_j g_k).
Generalised eigvec in mpmath, rounded to rationals, Rayleigh quotient Q(v)/||v||^2 evaluated in Arb => rigorous upper bound
for inf Q(f)/||f||^2 over the sector (Ritz/min-max). Usage: ritz_l2.py M L dps prec"""
import sys, json, time; sys.path.insert(0,'.'); A_=sys.argv[:]; sys.argv=sys.argv[:1]
from s46_planted import *
sys.argv=A_
M=int(sys.argv[1]); h=sys.argv[2]; mp.mp.dps=int(sys.argv[3]); prec=int(sys.argv[4])
t0=time.time(); A,T,ns=build_all(M,h,prec); ctx.prec=prec; ha=arb(h)
Gs=Gpolys(M)
Wz=arb_mat(M,M); Gr=arb_mat(M,M)
for i in range(M):
    for j in range(M):
        Wz[i,j]=A[i,j]-sum((2*arb(vonmangoldt(n)).log()/arb(n).sqrt())*T[n][i,j] for n in ns)
        if (i+j)%2==0:
            c=Gs[i,j].coeffs(); Gr[i,j]=ha*(arb(c[0].p)/c[0].q)
row={"M":M,"L":h,"build_sec":round(time.time()-t0,1)}
for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
    W=tomp(Wz,I); G=tomp(Gr,I)
    Lc=mp.cholesky(G); Li=mp.inverse(Lc); C=Li*W*Li.T; ev,U=mp.eigsy(C); V=Li.T*U
    k=min(range(len(I)),key=lambda t:ev[t]); v=V[:,k]
    sc=max(abs(t) for t in v); xs=[fmpq(int(mp.nint(t/sc*mp.mpf(10)**80)),10**80) for t in v]; xa=[arb(q.p)/q.q for q in xs]
    num=sum(xa[a]*xa[b]*Wz[I[a],I[b]] for a in range(len(I)) for b in range(len(I)))
    den=sum(xa[a]*xa[b]*Gr[I[a],I[b]] for a in range(len(I)) for b in range(len(I)))
    rq=num/den
    row[nm+"_ritz_mp"]=mp.nstr(ev[k],8); row[nm+"_RQ_arb_upper"]=rq.str(8); row[nm+"_RQ_upper_float"]=float(rq.upper())
row["sec"]=round(time.time()-t0,1); print(json.dumps(row),flush=True)
