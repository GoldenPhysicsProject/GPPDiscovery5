# NUMERICAL (SDP, CLARABEL): P.17.5 dual problem for the coupled line cost.
# c_min = sup{ 2c2 G(1,0) + 2c3 G(0,1) : G positive definite on Z^2, G(0)=1, G(m,n)=0 if |m a + n b| >= 2L }  (Prop P.17.D)
# Finite-section relaxation: only the box-K Toeplitz section is required PSD  =>  value_K >= c_min, nonincreasing in K.
import numpy as np, cvxpy as cp, itertools
a=np.log(2); b=np.log(3); c2=a/np.sqrt(2); c3=b/np.sqrt(3); c0=0.5*(c3+np.sqrt(c3**2+4*c2**2))
print('# p17_strip_sdp.py -- NUMERICAL. c0* = %.9f, c2+c3 = %.9f. RH OPEN.'%(c0,c2+c3))
for twoL,Ks in [(np.log(3)+0.05,[2,3,4,5,6,7,8]),(np.log(3)+0.2,[2,3,4,5]),(np.log(4)-0.01,[2,3,4,5])]:
    for K in Ks:
        P=list(itertools.product(range(-K,K+1),repeat=2))
        keys={}
        def var(k):
            k=tuple(k); kk=k if k>(0,0) else tuple(-x for x in k)
            return kk
        D=[(m,n) for m in range(-2*K,2*K+1) for n in range(-2*K,2*K+1)]
        free=sorted({var(k) for k in D if k!=(0,0) and abs(k[0]*a+k[1]*b)<twoL})
        idx={k:i for i,k in enumerate(free)}
        g=cp.Variable(len(free))
        n=len(P); Ms=[]
        rows=[]
        T=[[None]*n for _ in range(n)]
        for i,p in enumerate(P):
            for j,q in enumerate(P):
                k=(p[0]-q[0],p[1]-q[1])
                if k==(0,0): T[i][j]=1.0
                else:
                    kk=var(k); T[i][j]=g[idx[kk]] if kk in idx else 0.0
        Tm=cp.bmat([[cp.reshape(cp.hstack([T[i][j]]) if not isinstance(T[i][j],float) else cp.Constant(np.array([T[i][j]])),(1,1),order='C') for j in range(n)] for i in range(n)])
        obj=2*c2*(g[idx[(1,0)]] if (1,0) in idx else 0)+2*c3*(g[idx[(0,1)]] if (0,1) in idx else 0)
        pr=cp.Problem(cp.Maximize(obj),[ (Tm+Tm.T)/2>>0 ])
        pr.solve(solver='CLARABEL')
        print('2L=%.4f K=%d visible lags=%d  relaxation value=%.6f  (minus c0*: %+.2e)'%(twoL,K,len(free),pr.value,pr.value-c0),flush=True)
