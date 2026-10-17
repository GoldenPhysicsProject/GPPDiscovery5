"""Grok S62 H62a: rigorous L2-normalised Ritz UPPER bounds for the zeta Weil form per parity sector on [-L,L],
basis g_j(u)=P_j(u/L)*(1-(u/L)^2)^w  (w=4 is the S14/nogo_attack basis; w=0 = plain Legendre, jumps at +-L allowed:
F=g*g~ is then Lipschitz and the S14 real-space explicit formula is valid; g is in the Weil form domain).
Same exact assembly as nogo_attack/ritz_l2.py (S14/S44 Arb), Gram_jk = L*G_jk(0) = <g_j,g_k>_{L2}.
Generalised eigvec (mpmath) -> rationals (80 digits) -> Arb Rayleigh quotient Q(v)/||v||^2  = rigorous upper bound
for inf Q/||f||^2 over the parity sector (min-max).  Usage: ritz_w.py M L dps prec w [sectors]"""
import sys, json, time; sys.path.insert(0,'.'); A_=sys.argv[:]; sys.argv=sys.argv[:1]
import s14_exact
M=int(A_[1]); h=A_[2]; dps=int(A_[3]); prec=int(A_[4]); w=int(A_[5]); secs=A_[6] if len(A_)>6 else "even,odd"
Ww=s14_exact.fmpq_poly([1,0,-1])**w
s14_exact.p=lambda j: s14_exact.legendre(j)*Ww
from s46_planted import *
mp.mp.dps=dps
t0=time.time(); A,T,ns=build_all(M,h,prec); ctx.prec=prec; ha=arb(h)
Gs=Gpolys(M)
Wz=arb_mat(M,M); Gr=arb_mat(M,M)
for i in range(M):
    for j in range(M):
        Wz[i,j]=A[i,j]-sum((2*arb(vonmangoldt(n)).log()/arb(n).sqrt())*T[n][i,j] for n in ns)
        if (i+j)%2==0:
            c=Gs[i,j].coeffs(); Gr[i,j]=ha*(arb(c[0].p)/c[0].q)
row={"M":M,"L":h,"w":w,"dps":dps,"prec":prec,"build_sec":round(time.time()-t0,1),
     "maxrad_W":max(float(Wz[i,j].rad()) for i in range(M) for j in range(M))}
vecs={}
for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
    if nm not in secs: continue
    W=tomp(Wz,I); Gm=tomp(Gr,I)
    Lc=mp.cholesky(Gm); Li=mp.inverse(Lc); C=Li*W*Li.T; C=(C+C.T)/2; ev,U=mp.eigsy(C); V=Li.T*U
    o=sorted(range(len(I)),key=lambda t:ev[t]); k=o[0]; v=V[:,k]
    sc=max(abs(t) for t in v); xs=[fmpq(int(mp.nint(t/sc*mp.mpf(10)**80)),10**80) for t in v]; xa=[arb(q.p)/q.q for q in xs]
    num=sum(xa[a]*xa[b]*Wz[I[a],I[b]] for a in range(len(I)) for b in range(len(I)))
    den=sum(xa[a]*xa[b]*Gr[I[a],I[b]] for a in range(len(I)) for b in range(len(I)))
    rq=num/den
    row[nm+"_ritz_mp"]=mp.nstr(ev[k],10); row[nm+"_ritz2_mp"]=mp.nstr(ev[o[1]],10)
    row[nm+"_RQ_arb"]=rq.str(10); row[nm+"_RQ_upper"]=float(rq.upper())
    vecs[nm]=[str(q) for q in xs]
row["sec"]=round(time.time()-t0,1); print(json.dumps(row),flush=True)
json.dump({"row":row,"vecs":vecs},open(f"vec_M{M}_L{h}_w{w}.json","w"))
