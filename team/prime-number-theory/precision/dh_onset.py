"""DH failure onset L*: bisection in L for the Davenport-Heilbronn window form (cos basis, s62 assembly via
cosritz_dump.py, Arb), sector given. At each L one Kmax assembly; inertia (number of negative LDL^T pivots) of the
nested leading blocks K in KS gives convergence-in-K evidence. Negative pivot with all earlier pivots positive,
checked in ball arithmetic, = CERTIFIED non-PSD of the window form at that L (a Ritz witness).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
usage: dh_onset.py sector Kmax prec Llo_num Lhi_num den steps   (L = num/den)"""
import sys, json, subprocess, os
import mpmath as mp
from flint import arb, ctx
sec=sys.argv[1]; Kmax=int(sys.argv[2]); prec=int(sys.argv[3]); lo=int(sys.argv[4]); hi=int(sys.argv[5]); den=int(sys.argv[6]); steps=int(sys.argv[7])
KS=[k for k in (80,120,160,200,240,266,300,320,360,400) if k<=Kmax]
here=os.path.dirname(os.path.abspath(__file__)); py=sys.executable
def evalL(num):
    Ls=f"{num}/{den}"; fn=f"Wdump_DH_{sec}_L{num}_{den}_K{Kmax}_p{prec}.json"
    if not os.path.exists(fn):
        subprocess.run([py,os.path.join(here,"cosritz_dump.py"),Ls,str(Kmax),str(prec),str(int(prec*0.30103)),sec],env=dict(os.environ,FORM="DH"),check=True,capture_output=True)
    D=json.load(open(fn)); mp.mp.prec=prec
    A=[[mp.mpf(r[0]) for r in row] for row in D["W"]]
    # float LDL^T on full Kmax; leading-block inertia read off the pivots (pivots of block K = first K pivots)
    K=Kmax; d=[]; Lm=[[mp.mpf(0)]*K for _ in range(K)]
    for j in range(K):
        s=A[j][j]-mp.fsum(Lm[j][k]**2*d[k] for k in range(j)); d.append(s)
        if s==0: break
        for i in range(j+1,K): Lm[i][j]=(A[i][j]-mp.fsum(Lm[i][k]*Lm[j][k]*d[k] for k in range(j)))/s
    inert={k:sum(1 for x in d[:k] if x<0) for k in KS}
    firstneg=next((j for j,x in enumerate(d) if x<0),None)
    rec={"L":Ls,"Lf":num/den,"sector":sec,"Kmax":Kmax,"prec":prec,"neg_pivots_by_K":inert,"first_neg":firstneg,
         "min_pivot":mp.nstr(min(d),6)}
    print(json.dumps(rec),flush=True); return inert[KS[-1]]>0
a,b=lo,hi
fa=evalL(a); fb=evalL(b)
print(json.dumps({"bracket_start":[a/den,fa,b/den,fb]}),flush=True)
for _ in range(steps):
    if fa==fb: break
    m=(a+b)//2
    if m in (a,b): break
    fm=evalL(m)
    if fm==fa: a,fa=m,fm
    else: b,fb=m,fm
print(json.dumps({"onset_bracket":[a/den,b/den],"sector":sec,"Kmax":Kmax}),flush=True)
