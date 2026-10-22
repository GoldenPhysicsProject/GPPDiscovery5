"""Grok S80 H80b: error-free sliced (Ozaki-type) Gram product G = V diag(w) V^T for 320-bit fixed-point V, using only
float64 dgemm with exact integer slices (each slice b bits, k terms: need 2b + log2 k <= 53). Checks exactness against
Python big integers and times dgemm throughput. NUMERICAL timing; exactness is by construction (integer dgemm < 2^53)."""
import numpy as np, time, random, math, sys
P=320; k=64; b=(53-math.ceil(math.log2(k)))//2   # slice bits
S=math.ceil(P/b)
def slices(M):   # M: object array of ints in [0,2^P) (signs handled separately)
    out=[]; 
    for s in range(S):
        sh=P-(s+1)*b
        out.append(np.vectorize(lambda x: float((x>>max(sh,0)) & ((1<<b)-1) if sh>=0 else (x & ((1<<(b+sh))-1))<< (-sh) ))(M).astype(np.float64))
    return out
def test_exact(n=40):
    random.seed(1)
    V=np.array([[random.getrandbits(P) for _ in range(k)] for _ in range(n)],dtype=object)
    Vs=slices(V)
    # recombine exactly: sum over i,j of 2^{2P-(i+1)b-(j+1)b} * (Vs_i Vs_j^T)
    G=np.zeros((n,n),dtype=object)
    for i in range(S):
        for j in range(S):
            C=(Vs[i]@Vs[j].T); assert np.abs(C).max()<2**53
            sh=2*P-(i+1)*b-(j+1)*b
            Ci=C.astype(np.int64).astype(object)
            G+= (Ci<<sh) if sh>=0 else (Ci>>(-sh))
    Gex=V.dot(V.T)
    d=np.abs((Gex-G)).max(); return d, S, b
d,S_,b_=test_exact()
print('slice bits',b_,'slices',S_,'max|exact-sliced| =',d,'(nonzero only from final-slice truncation if P % b)')
Np=int(sys.argv[1]) if len(sys.argv)>1 else 3000
A=np.random.randint(0,2**b,(Np,k)).astype(np.float64)
t=time.time(); 
for _ in range(5): C=A@A.T
t=(time.time()-t)/5
flop=2*Np*Np*k; print('Np=%d k=%d one slice dgemm %.3fs -> %.1f GFLOP/s; pairs per cell job = %d -> %.1fs per job (vs Arb 320-bit 14.6s at N2732)'%(Np,k,t,flop/t/1e9,S*(S+1)//2,t*S*(S+1)//2))
