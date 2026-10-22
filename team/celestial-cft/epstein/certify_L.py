"""CERTIFIED zero counts for Z_Q (D=-4,-20,-84; genus-only class groups) via rigorous argument principle.
Each boundary segment S is enclosed in a complex ball; Z_Q is evaluated on that ball with Arb (zq_L.py).
If the enclosure excludes 0 and lies in a half-plane, arg change on S = Arg(Z(end)/Z(start)) (rigorous);
otherwise S is bisected.  Total/2pi must enclose exactly one integer.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, time
sys.path.insert(0,'.')
from zq_L import ZQ
from flint import acb, arb, ctx
from multiprocessing import Pool
ctx.prec=80

def ev(Q,p): return ZQ(Q,acb(p.real,p.imag))
def evball(Q,p,q):
    m=(p+q)/2
    return ZQ(Q,acb(arb(m.real,abs((q-p).real)/2+1e-30),arb(m.imag,abs((q-p).imag)/2+1e-30)))
def ok(W): return (W.real>0) or (W.real<0) or (W.imag>0) or (W.imag<0)

def seg(Q,p,q,fp,fq,d,stats):
    W=evball(Q,p,q); stats[0]+=1
    if ok(W):
        return (fq/fp).arg()
    if d>40: raise RuntimeError('segment refinement failed near %s'%p)
    m=(p+q)/2; fm=ev(Q,m)
    return seg(Q,p,m,fp,fm,d+1,stats)+seg(Q,m,q,fm,fq,d+1,stats)

def count(Q,x0,x1,y0,y1,n0=40):
    c=[complex(x0,y0),complex(x1,y0),complex(x1,y1),complex(x0,y1)]
    pts=[]
    for i in range(4):
        a,b=c[i],c[(i+1)%4]
        pts+=[a+(b-a)*k/n0 for k in range(n0)]
    pts.append(pts[0])
    vals=[ev(Q,p) for p in pts]
    tot=arb(0); stats=[0]
    for i in range(len(pts)-1):
        tot+=seg(Q,pts[i],pts[i+1],vals[i],vals[i+1],0,stats)
    k=tot/(2*arb.pi()); n=int(float(k.mid())+0.5)
    assert (k-n).abs_upper()<0.5, k
    return n, str(k), stats[0]

def job(a):
    Q,box=a; t=time.time()
    n,k,ns=count(Q,*box)
    return dict(Q=Q,box=box,count=n,count_ball=k,ball_evals=ns,sec=round(time.time()-t,1))

if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='big':   # region [0.52,3] x [1,T]
        T=float(sys.argv[2]); forms=[(1,0,1),(1,0,5),(2,2,3),(1,0,21),(3,0,7),(2,2,11),(5,4,5)]
        tasks=[]
        for Q in forms:
            y=1.0
            while y<T: tasks.append((Q,(0.52,3.0,y,min(y+10,T)))); y+=10
        with Pool(int(sys.argv[3]) if len(sys.argv)>3 else 4) as p: res=p.map(job,tasks,chunksize=1)
        json.dump(res,open(f'certified_counts_T{int(T)}.json','w'),indent=1)
        for r in res: print(r)
    elif mode=='zeros':  # small boxes around located zeros from file
        zs=json.load(open(sys.argv[2])); r=float(sys.argv[3]) if len(sys.argv)>3 else 1e-8
        tasks=[(tuple(z['Q']),(z['x']-r,z['x']+r,z['y']-r,z['y']+r)) for z in zs]
        with Pool(4) as p: res=p.map(job,tasks,chunksize=1)
        json.dump(res,open('certified_zero_boxes.json','w'),indent=1)
        for r_ in res: print(r_)
