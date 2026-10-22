"""lambda_max of the prime operator C_p on L^2(-L,L) (half-width L), Galerkin in the ON cos/sin basis, float64 (NUMERICAL).
(C_p f)(x) = sum_{n<e^{2L}} Lambda(n) n^{-1/2} [f(x+log n) + f(x-log n)], f extended by 0. Feeds the Astra/Mori strip
bound Omega > 2 pi e^{s+m}. Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
usage: cp_norm.py L K1,K2,.."""
import sys, numpy as np, math, json
L=float(eval(sys.argv[1])); Ks=[int(k) for k in sys.argv[2].split(',')]
def mang(n):
    for p in range(2,n+1):
        if n%p==0:
            m=n
            while m%p==0: m//=p
            return math.log(p) if m==1 else 0.0
ns=[n for n in range(2,int(math.exp(2*L))+1) if mang(n)>0 and math.log(n)<2*L]
def basis(x,K,par):
    if par=='even': return np.array([np.cos((2*m+1)*np.pi*x/(2*L)) for m in range(K)]).T/np.sqrt(L)
    return np.array([np.sin((m+1)*np.pi*x/L) for m in range(K)]).T/np.sqrt(L)
out={"L":L,"ns":ns}
for K in Ks:
    nq=max(400,int(3*K)); gx,gw=np.polynomial.legendre.leggauss(nq)
    for par in ('even','odd'):
        M=np.zeros((K,K))
        for n in ns:
            a=math.log(n); w=mang(n)/math.sqrt(n)
            lo,hi=-L,L-a; x=(hi-lo)/2*gx+(hi+lo)/2; ww=gw*(hi-lo)/2
            A=basis(x,K,par); B=basis(x+a,K,par)
            S=(A*ww[:,None]).T@B          # <phi_m, f(.+a)> component: (S_a phi_n)(x)=phi_n(x+a) on x in [-L,L-a]
            M+=w*(S+S.T)                  # S_{-a} is the adjoint
        ev=np.linalg.eigvalsh(M); out[f"{par}_K{K}"]=[float(ev[-1]),float(ev[0])]
print(json.dumps(out))
