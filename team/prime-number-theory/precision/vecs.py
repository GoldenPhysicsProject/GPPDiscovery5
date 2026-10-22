"""Smallest-|eigenvalue| eigenvector (the soft / onset mode) of a dumped K-section by inverse iteration on an LDL^T
factorisation (works for the indefinite DH section too). Unnormalised cos basis f = sum c_m cos(a_m u), a_m=(2m+1)pi/(2L);
output coefficients are rescaled so that ||f||_{L^2(-L,L)} = 1. Also endpoint derivatives f'(L), f'''(L) (decay of f^).
NUMERICAL. Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
usage: vecs.py dump.json K out.json [prec]"""
import sys, json, mpmath as mp
D=json.load(open(sys.argv[1])); K=int(sys.argv[2]); out=sys.argv[3]
mp.mp.prec=int(sys.argv[4]) if len(sys.argv)>4 else D["prec"]
p,q=D["L"].split('/'); L=mp.mpf(p)/int(q)
A=[[mp.mpf(D["W"][i][j][0]) for j in range(K)] for i in range(K)]
d=[]; Lm=[[mp.mpf(0)]*K for _ in range(K)]
for j in range(K):
    s=A[j][j]-mp.fsum(Lm[j][k]**2*d[k] for k in range(j)); d.append(s)
    for i in range(j+1,K): Lm[i][j]=(A[i][j]-mp.fsum(Lm[i][k]*Lm[j][k]*d[k] for k in range(j)))/s
def solve(b):
    y=list(b)
    for i in range(K): y[i]=y[i]-mp.fsum(Lm[i][k]*y[k] for k in range(i))
    y=[y[i]/d[i] for i in range(K)]
    for i in reversed(range(K)): y[i]=y[i]-mp.fsum(Lm[k][i]*y[k] for k in range(i+1,K))
    return y
v=[mp.mpf(1)/(i+1) for i in range(K)]
for it in range(30):
    w=solve(v); n=mp.sqrt(mp.fsum(x*x for x in w)); v=[x/n for x in w]
Av=[mp.fsum(A[i][j]*v[j] for j in range(K)) for i in range(K)]
lamW=mp.fsum(v[i]*Av[i] for i in range(K)); resid=mp.sqrt(mp.fsum((Av[i]-lamW*v[i])**2 for i in range(K)))
a=[(2*m+1)*mp.pi/(2*L) for m in range(K)]
c=[x/mp.sqrt(L) for x in v]   # ||f||^2 = L*sum v^2 = 1
d1=mp.fsum(-c[m]*a[m]*mp.sin(a[m]*L) for m in range(K)); d3=mp.fsum(c[m]*a[m]**3*mp.sin(a[m]*L) for m in range(K))
rec={"form":D["form"],"L":D["L"],"K":K,"Q_over_norm2":mp.nstr(lamW/L,12),"eig_W":mp.nstr(lamW,12),"residual_W":mp.nstr(resid,4),
     "fp_L":mp.nstr(d1,8),"fppp_L":mp.nstr(d3,8),"coef":[mp.nstr(x,int(mp.mp.prec*0.3)) for x in c]}
json.dump(rec,open(out,"w")); print(json.dumps({k:v for k,v in rec.items() if k!="coef"}))
