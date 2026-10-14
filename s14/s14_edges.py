"""S14 edges: Ritz window J_M(h) for A - cT on even/odd sectors, high precision (mpmath dps=120 from Arb entries), M up to 40.
Then certify: take null vector x of A - c_edge T (rounded to rationals), compute c* = <x,Ax>/<x,Tx> in Arb.
 even: <x,Tx> > 0  =>  A - cT not PSD for all c > c*  (rigorous upper bound on true even edge).
 odd : <x,Tx> < 0  =>  A - cT not PSD for all c < c*  (rigorous lower bound on true odd edge)."""
import sys, json; sys.path.insert(0,'.')
from s14_exact import *; import mpmath as mp
mp.mp.dps=120
def tomp(Ma,I): return mp.matrix([[mp.mpf(Ma[i,j].mid().str(150,radius=False)) for j in I] for i in I])
def psd(X):
    try: mp.cholesky(X); return True
    except Exception: return False
M=int(sys.argv[1]); hs=sys.argv[2].split(',')
cz=mp.sqrt(2)*mp.log(2); czA=arb(2).sqrt()*arb(2).log(); out=[]
for hstr in hs:
    A,T=build(M,hstr,int(sys.argv[3]) if len(sys.argv)>3 else 768)
    maxrad=max(float(A[i,j].rad()) for i in range(M) for j in range(M))
    row={"M":M,"h":hstr,"maxrad":maxrad}
    for nm,I,sgn in (("even",list(range(0,M,2)),1),("odd",list(range(1,M,2)),-1)):
        a,t=tomp(A,I),tomp(T,I)
        if not psd(a-cz*t): row[nm]="c_zeta NOT in Ritz window"; continue
        lo,hi=cz,cz+sgn*mp.mpf('0.05')
        for _ in range(110):
            m=(lo+hi)/2
            if psd(a-m*t): lo=m
            else: hi=m
        ce=lo
        ev,V=mp.eigsy(a-ce*t); x=V[:,0]
        xs=[fmpq(int(mp.nint(x[i]*mp.mpf(10)**60)),10**60) for i in range(len(I))]
        av=sum(arb(xs[i].p)/xs[i].q*arb(xs[j].p)/xs[j].q*A[I[i],I[j]] for i in range(len(I)) for j in range(len(I)))
        tv=sum(arb(xs[i].p)/xs[i].q*arb(xs[j].p)/xs[j].q*T[I[i],I[j]] for i in range(len(I)) for j in range(len(I)))
        cs=av/tv; gap=cs-czA
        row[nm]={"ritz_edge":mp.nstr(ce,25),"edge_minus_czeta":mp.nstr(ce-cz,12),"cert_cstar":cs.str(25),"cert_cstar_minus_czeta":gap.str(10),
                 "tx_sign":"+" if tv>0 else ("-" if tv<0 else "?"),"W_zeta(x)_over_|tx|":gap.str(8)}
    print(json.dumps(row),flush=True); out.append(row)
