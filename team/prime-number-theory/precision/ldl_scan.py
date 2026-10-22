"""Precision scan for the LDL^T factorization of the zeta Weil window matrix (cos basis, even sector),
with Davenport-Heilbronn at the same L as the failing control.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

Input: Wdump_*.json from cosritz_dump.py (entries assembled in Arb at high precision, midpoint + radius).
For each basis size K (leading K x K block; the cos basis is nested so this IS the K-mode matrix) and each
working precision p (bits):
  (a) FLOAT: round the midpoints to p bits, LDL^T in p-bit floating point (mpmath). Record whether all
      pivots are > 0, the min pivot, and its relative error against the reference (dps of the dump).
  (b) BALL: entries as Arb balls at p bits with radius max(dump radius, 2^(8-p)*|entry|) (models assembly
      at p bits with ~8 bits of assembly loss), ball LDL^T at p bits. CERTIFIED PD of the finite K-matrix
      iff every pivot ball is > 0; CERTIFIED indefinite iff some pivot ball is < 0 with all earlier > 0.
NOTE: PD of the finite K-section is a necessary condition only; it is NOT a Weil-positivity certificate
at L (that needs a tail/truncation bound, e.g. the Mori route). RH is OPEN.
usage: ldl_scan.py dumpfile "K1,K2,.." "p1,p2,.." [ball]
"""
import sys, json, time, math
import mpmath as mp
from flint import arb, ctx
D=json.load(open(sys.argv[1])); Ks=[int(x) for x in sys.argv[2].split(',')]; ps=[int(x) for x in sys.argv[3].split(',')]
doball=len(sys.argv)>4
Kmax=D["K"]; Wm=[[r[0] for r in row] for row in D["W"]]; Wr=[[r[1] for r in row] for row in D["W"]]

def ldl_float(K,p):
    mp.mp.prec=p
    A=[[mp.mpf(Wm[i][j]) for j in range(K)] for i in range(K)]
    d=[]; neg=None
    Lm=[[mp.mpf(0)]*K for _ in range(K)]
    for j in range(K):
        s=A[j][j]-mp.fsum(Lm[j][k]**2*d[k] for k in range(j))
        d.append(s)
        if s<=0 and neg is None: neg=j
        if s==0: break
        for i in range(j+1,K):
            Lm[i][j]=(A[i][j]-mp.fsum(Lm[i][k]*Lm[j][k]*d[k] for k in range(j)))/s
    return d,neg

def ldl_ball(K,p):
    ctx.prec=p
    A=[[arb(Wm[i][j]) for j in range(K)] for i in range(K)]
    floor=arb(2)**(8-p)
    for i in range(K):
        for j in range(K):
            r=max(arb(Wr[i][j]), floor*abs(A[i][j]).upper())
            A[i][j]=arb(A[i][j].mid(), r.upper())
    d=[]; Lm=[[None]*K for _ in range(K)]
    for j in range(K):
        s=A[j][j]
        for k in range(j): s-=Lm[j][k]*Lm[j][k]*d[k]
        d.append(s)
        if not s>0:
            return d, ("NEG" if s<0 else "UNRESOLVED"), j
        for i in range(j+1,K):
            t=A[i][j]
            for k in range(j): t-=Lm[i][k]*Lm[j][k]*d[k]
            Lm[i][j]=t/s
    return d,"PD",K

refp=int(D.get("prec",960))
for K in Ks:
    t=time.time(); dref,negref=ldl_float(K,refp); mp.mp.prec=refp
    jmin=min(range(len(dref)),key=lambda j: dref[j]) if negref is None else negref
    print(json.dumps({"form":D["form"],"L":D["L"],"K":K,"ref_prec":refp,"ref_first_nonpos":negref,
        "ref_min_pivot":mp.nstr(dref[jmin],8),"at":jmin,"ref_last_pivot":mp.nstr(dref[-1],8),
        "max_dump_rad":max(max(r[:K]) for r in Wr[:K]),"sec":round(time.time()-t,1)}),flush=True)
    for p in ps:
        t=time.time(); d,neg=ldl_float(K,p); mp.mp.prec=refp
        rec={"K":K,"p":p,"float_first_nonpos":neg}
        if neg is None or negref is not None:
            j=jmin if negref is None else negref
            if j<len(d): rec["pivot"]=mp.nstr(mp.mpf(d[j]),6); rec["relerr"]=mp.nstr(abs((mp.mpf(d[j])-dref[j])/dref[j]),3)
        rec["fsec"]=round(time.time()-t,1)
        if doball:
            t=time.time(); db,st,j=ldl_ball(K,p); rec.update({"ball":st,"ball_stop":j,"ball_piv":db[-1].str(5) if db else None,"bsec":round(time.time()-t,1)})
        print(json.dumps(rec),flush=True)
