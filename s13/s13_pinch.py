"""S13b: width of the admissible coupling window for prime 2 (zeta pole+arch), bisection. NUMERICAL float Ritz M=14."""
import numpy as np, sys
from scipy.linalg import eigh
exec(open("s13/s13_coupling.py").read().split("for th in")[0])
def edge(a,t,g,c0,dirn):
    lo,hi=c0,c0+dirn*3
    f=lambda c: eigh(a-c*t,g,eigvals_only=True)[0]>=0
    if not f(c0): return None
    for _ in range(50):
        m=(lo+hi)/2; (lo,hi)=(m,hi) if f(m) else (lo,m)
    return lo
c=np.sqrt(2)*L2
for h in (0.50,0.52,0.53,0.54,0.545,0.549):
    Q,G=ns["form"](h,"Z",0.0); T=T2(h); A=Q+c*T; r=[]
    for I in (list(range(0,M,2)),list(range(1,M,2))):
        a,t,g=A[np.ix_(I,I)],T[np.ix_(I,I)],G[np.ix_(I,I)]; r+= [edge(a,t,g,c,-1),edge(a,t,g,c,+1)]
    lo=max(r[0],r[2]); hi=min(r[1],r[3])
    print(f"h={h:.3f} even[{r[0]:.6f},{r[1]:.6f}] odd[{r[2]:.6f},{r[3]:.6f}]  joint[{lo:.6f},{hi:.6f}] width={hi-lo:.2e}  c_zeta={c:.6f}  c_F1/4={(2**.25+2**-.25)*L2/np.sqrt(2):.6f}",flush=True)
