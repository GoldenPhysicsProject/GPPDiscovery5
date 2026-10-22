"""NUMERICAL search for zeros of Z_{x^2+5y^2} with Re s > 1 (lead S100 item 2); certification is done separately.
For Re s > 1, Z = L(s,chi_-4) L(s,chi_5) (1 + R(s)) with R = zeta L(s,chi_-20) / (L(s,chi_-4) L(s,chi_5)), and the L-factor has no zeros.
So zeros with Re s > sigma0 are zeros of 1 + R. Their number in [sigma0, 3] x [T1, T2] is the winding of 1 + R along the boundary.
On Re s = 3, |R - R(3)| is small and 1 + R stays away from 0.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, cmath, math
sys.path.insert(0,'.')
from zq_L import Ld
from flint import acb, ctx
from multiprocessing import Pool
ctx.prec=53
def G(s):
    s=acb(s.real,s.imag)
    v=1+ s.zeta()*Ld(-20,s)/(Ld(-4,s)*Ld(5,s))
    return complex(float(v.real.mid()),float(v.imag.mid()))
def wind(pts):
    vals=[G(p) for p in pts]; w=0.0
    for a,b in zip(vals,vals[1:]): w+=cmath.phase(b/a)
    return w/(2*math.pi), vals
def box(args):
    s0,T1,T2,h=args
    n=int((T2-T1)/h)+1
    right=[complex(s0,T1+k*(T2-T1)/n) for k in range(n+1)]
    top=[complex(s0+k*(3-s0)/200,T2) for k in range(201)]
    left=[complex(3,T2-k*(T2-T1)/n*10) for k in range(n//10+1)]
    bot=[complex(3-k*(3-s0)/200,T1) for k in range(201)]
    w=0
    for seg in (right,top,left,bot): w+=wind(seg)[0]
    return (T1,T2,round(-w,3))   # boundary traversed clockwise (up the left edge s0), hence the sign
if __name__=='__main__':
    s0=float(sys.argv[1]); T=float(sys.argv[2]); h=float(sys.argv[3]); W=int(sys.argv[4])
    tasks=[(s0,t,t+10,h) for t in range(1,int(T),10)]
    out=[]
    with Pool(W) as p:
        for r in p.imap(box,tasks):
            if abs(r[2])>0.5: print(r,flush=True)
            out.append(r)
    json.dump(out,open(f'search_re_gt1_s{s0}_T{int(T)}.json','w'))
    print('done',flush=True)
