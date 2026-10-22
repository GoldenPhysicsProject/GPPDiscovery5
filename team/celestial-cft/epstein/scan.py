"""Locate zeros of Epstein zeta Z_Q(s) with Re s > SIG0 by the argument principle (float-level,
adaptive sampling = NUMERICAL), then Newton-refine.  Certification is done separately in certify.py.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, math, cmath
sys.path.insert(0,'.')
from theta import Lambda
from flint import acb, ctx
from multiprocessing import Pool

SIG0, SIG1 = 0.52, 3.0
FORMS = {
 -4:  [(1,0,1)],
 -20: [(1,0,5),(2,2,3)],
 -23: [(1,1,6),(2,1,3)],
 -56: [(1,0,14),(2,0,7),(3,2,5)],
 -84: [(1,0,21),(3,0,7),(2,2,11),(5,4,5)],
}

def F(Q,s):
    ctx.prec=53
    v=Lambda(Q,acb(s.real,s.imag),bits=26)
    return complex(float(v.real.mid()),float(v.imag.mid()))

def edge_winding(Q,a,b,h0=0.05,depth=0):
    # accumulate arg change along segment a->b with adaptive bisection
    n=max(2,int(abs(b-a)/h0))
    pts=[a+(b-a)*k/n for k in range(n+1)]
    vals=[F(Q,p) for p in pts]
    tot=0.0
    for i in range(n):
        tot+=seg(Q,pts[i],pts[i+1],vals[i],vals[i+1],0)
    return tot

def seg(Q,p,q,fp,fq,d):
    dth=cmath.phase(fq/fp)
    if (abs(dth)>0.4 or abs(math.log(abs(fq/fp)))>1.0) and d<14:
        m=(p+q)/2; fm=F(Q,m)
        return seg(Q,p,m,fp,fm,d+1)+seg(Q,m,q,fm,fq,d+1)
    return dth

def winding(Q,x0,x1,y0,y1):
    c=[complex(x0,y0),complex(x1,y0),complex(x1,y1),complex(x0,y1)]
    tot=sum(edge_winding(Q,c[i],c[(i+1)%4]) for i in range(4))
    return round(tot/(2*math.pi)), tot/(2*math.pi)

def locate(Q,x0,x1,y0,y1,k,out):
    if k==0: return
    if (x1-x0)<0.02 and (y1-y0)<0.02:
        out.append(((x0+x1)/2,(y0+y1)/2,k)); return
    xm=(x0+x1)/2; ym=(y0+y1)/2
    for (a,b,c,d) in [(x0,xm,y0,ym),(xm,x1,y0,ym),(x0,xm,ym,y1),(xm,x1,ym,y1)]:
        kk,_=winding(Q,a,b,c,d)
        locate(Q,a,b,c,d,kk,out)

def newton(Q,z):
    import mpmath
    f=lambda s: F(Q,complex(s))
    for _ in range(40):
        h=1e-6
        fz=f(z); d=(f(z+h)-f(z-h))/(2*h)
        dz=fz/d; z=z-dz
        if abs(dz)<1e-13: break
    return z

def job(args):
    Q,y0,y1=args
    k,raw=winding(Q,SIG0,SIG1,y0,y1)
    found=[]
    if k>0:
        cells=[]; locate(Q,SIG0,SIG1,y0,y1,k,cells)
        for (x,y,m) in cells:
            z=newton(Q,complex(x,y)); found.append((z.real,z.imag,m))
    return dict(Q=Q,y0=y0,y1=y1,count=k,raw=raw,zeros=found)

if __name__=='__main__':
    T=float(sys.argv[1]) if len(sys.argv)>1 else 100
    tasks=[]
    for D,fs in FORMS.items():
        for Q in fs:
            y=1.0
            while y<T:
                tasks.append((Q,y,min(y+10,T))); y+=10
    with Pool(6) as p:
        res=p.map(job,tasks,chunksize=1)
    json.dump(res,open(f'scan_T{int(T)}.json','w'),indent=1)
    for r in res:
        if r['count'] or abs(r['raw']-round(r['raw']))>0.05:
            print(r['Q'],r['y0'],r['y1'],r['count'],round(r['raw'],3),r['zeros'])
    print('done')
