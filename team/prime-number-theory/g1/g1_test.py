"""Numerical test of Conjecture G1 (integer_regularity.md §4). NUMERICAL.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
System: rational primes with a random subset P_S deleted, P(p in P_S) = p^(alpha-1) independently (p >= p0).
Integers = n <= X coprime to P_S. Error E(x) = C(x) - D x, with D fitted by least squares (this absorbs the density,
including the contribution of deleted primes > X). Reports M_k = max_{2^k <= x < 2^{k+1}} |E(x)| and the fitted growth exponent.
Predictions: BDR hyperbola bound 2a/(a+2); Conjecture G1 a/2 (+eps).
usage: g1_test.py alpha X seed [p0]"""
import sys, numpy as np, json, math, time
a=float(sys.argv[1]); X=int(float(sys.argv[2])); seed=int(sys.argv[3]); p0=int(sys.argv[4]) if len(sys.argv)>4 else 2
t0=time.time(); rng=np.random.default_rng(seed)
s=np.ones(X//2+1,dtype=bool); s[0]=False          # odd sieve: index i <-> 2i+1
for i in range(1,int(math.isqrt(X))//2+1):
    if s[i]: p=2*i+1; s[p*p//2::p]=False
primes=np.concatenate(([2],2*np.nonzero(s)[0]+1)); primes=primes[primes<=X]; del s
pr=np.where(primes>=p0, primes.astype(float)**(a-1), 0.0)
S=primes[rng.random(len(primes))<pr]
keep=np.ones(X+1,dtype=bool); keep[0]=False
for p in S: keep[p::p]=False
# pass 1: sample C(x) to fit D
B=10**7; xs=[]; cs=[]; cnt=0
for st in range(0,X+1,B):
    blk=keep[st:st+B]; c=np.cumsum(blk,dtype=np.int64)+cnt
    idx=np.arange(0,len(blk),1000); xs.append(st+idx); cs.append(c[idx]); cnt=int(c[-1])
xs=np.concatenate(xs).astype(float); cs=np.concatenate(cs).astype(float)
m=xs>=1e3; D=float(np.sum(xs[m]*cs[m])/np.sum(xs[m]**2))
# pass 2: exact dyadic maxima
K=int(math.log2(X)); M=np.zeros(K+1); cnt=0
for st in range(0,X+1,B):
    blk=keep[st:st+B]; c=np.cumsum(blk,dtype=np.int64)+cnt; cnt=int(c[-1])
    x=np.arange(st,st+len(blk),dtype=float); E=np.abs(c-D*x)
    k=np.floor(np.log2(np.maximum(x,1))).astype(int)
    np.maximum.at(M,k,E)
ks=np.arange(K+1); sel=(ks>=12)&(ks<=K-1)
slope=np.polyfit(ks[sel]*math.log(2),np.log(M[sel]),1)[0]
selh=(ks>=K-8)&(ks<=K-1); slope_hi=np.polyfit(ks[selh]*math.log(2),np.log(M[selh]),1)[0]
print(json.dumps({"alpha":a,"X":X,"seed":seed,"p0":p0,"n_deleted":int(len(S)),"D":D,"slope_k12_up":round(float(slope),4),
  "slope_last8":round(float(slope_hi),4),"pred_BDR":round(2*a/(a+2),4),"pred_G1":a/2,"M":[float(v) for v in M],"sec":round(time.time()-t0,1)}))
