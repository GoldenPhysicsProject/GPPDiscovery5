"""Rigorous (Arb) positive-definiteness certificate for the K-mode zeta Weil window matrix from a p-bit
floating LDL^T, and scan of the smallest p for which it succeeds. DH at the same L is the failing control.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

Method (congruence + residual, Rump-style): A_ball = entries from the Arb assembly, radius inflated to
max(assembly radius, 2^(8-p)|a_ij|) so that it models an assembly done at p bits. Float LDL^T of mid(A) at p bits
gives exact binary L~ (unit lower), d~. Then A = L~ (D~ + G) L~^T with G = L~^{-1} E L~^{-T}, E = A - L~ D~ L~^T.
If all d~ > 0 and || D~^{-1/2} G D~^{-1/2} ||_F < 1 (computed in Arb), every matrix in A_ball is PD: CERTIFIED PD of
the finite K-section. (Necessary, not sufficient, for Weil positivity at L: no tail bound here. RH is OPEN.)
Also reports lambda_min of the K-section (inverse iteration, reference precision): NUMERICAL.
usage: ldl_certify.py dump.json "K1,.." "p1,.."  """
import sys, json, time
import mpmath as mp
from flint import arb, arb_mat, ctx, fmpz
def A2(x):
    sg,m,e,bc=mp.mpf(x)._mpf_
    if not m: return arb(0)
    v=arb(fmpz(int(m)))*arb(2)**int(e)
    return -v if sg else v
D=json.load(open(sys.argv[1])); Ks=[int(x) for x in sys.argv[2].split(',')]; ps=[int(x) for x in sys.argv[3].split(',')]
Wm=[[r[0] for r in row] for row in D["W"]]; Wr=[[r[1] for r in row] for row in D["W"]]; refp=D["prec"]
def fldl(K,p):
    mp.mp.prec=p
    A=[[mp.mpf(Wm[i][j]) for j in range(K)] for i in range(K)]
    d=[]; Lm=[[mp.mpf(0)]*K for _ in range(K)]
    for j in range(K):
        s=A[j][j]-mp.fsum(Lm[j][k]**2*d[k] for k in range(j)); d.append(s)
        if s<=0: return d,Lm,j
        for i in range(j+1,K): Lm[i][j]=(A[i][j]-mp.fsum(Lm[i][k]*Lm[j][k]*d[k] for k in range(j)))/s
    return d,Lm,None
def lam_min(K):
    d,Lm,bad=fldl(K,refp)
    if bad is not None: return None
    mp.mp.prec=refp
    def solve(b):
        y=list(b)
        for i in range(K): y[i]=y[i]-mp.fsum(Lm[i][k]*y[k] for k in range(i))
        y=[y[i]/d[i] for i in range(K)]
        for i in reversed(range(K)): y[i]=y[i]-mp.fsum(Lm[k][i]*y[k] for k in range(i+1,K))
        return y
    v=[mp.mpf(1)/(i+1) for i in range(K)]
    for it in range(40):
        w=solve(v); n=mp.sqrt(mp.fsum(x*x for x in w)); v=[x/n for x in w]
    Av=[mp.fsum(mp.mpf(Wm[i][j])*v[j] for j in range(K)) for i in range(K)]
    return mp.fsum(v[i]*Av[i] for i in range(K))
def certify(K,p):
    d,Lm,bad=fldl(K,p)
    if bad is not None: return {"cert":"FLOAT_LDL_FAILED","at":bad}
    ctx.prec=p+64
    floor=arb(2)**(8-p)
    A=arb_mat(K,K); Lt=arb_mat(K,K)
    for i in range(K):
        for j in range(K):
            a=arb(Wm[i][j]); r=max(arb(Wr[i][j]),floor*abs(a).upper()); A[i,j]=arb(a.mid(),r.upper())
            Lt[i,j]=arb(1) if i==j else (A2(Lm[i][j]) if j<i else arb(0))
    Dm=arb_mat(K,K)
    for i in range(K): Dm[i,i]=A2(d[i])
    E=A-Lt*Dm*Lt.transpose()
    # X = L^{-1} by forward substitution (unit lower, exact binary entries)
    X=arb_mat(K,K)
    for i in range(K):
        X[i,i]=arb(1)
        for j in range(i):
            s=arb(0)
            for k in range(j,i): s-=Lt[i,k]*X[k,j]
            X[i,j]=s
    S=arb_mat(K,K)
    for i in range(K): S[i,i]=1/A2(d[i]).sqrt()
    Y=S*X; G=Y*E*Y.transpose()
    fro=arb(0)
    for i in range(K):
        for j in range(K): fro+=abs(G[i,j]).upper()**2
    fro=fro.upper().sqrt().upper()
    ok=bool(fro<1)
    return {"cert":"PD" if ok else "NOT_CERTIFIED","fro_bound":fro.str(5),"min_pivot":mp.nstr(min(d),5)}
for K in Ks:
    t=time.time(); lm=lam_min(K)
    print(json.dumps({"form":D["form"],"L":D["L"],"K":K,"lambda_min_ref":mp.nstr(lm,8) if lm is not None else "LDL fails at ref (not PD)","sec":round(time.time()-t,1)}),flush=True)
    for p in ps:
        t=time.time(); r=certify(K,p); r.update({"K":K,"p":p,"sec":round(time.time()-t,1)}); print(json.dumps(r),flush=True)
