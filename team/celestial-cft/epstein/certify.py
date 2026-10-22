"""Rigorous (ball arithmetic, python-flint/Arb) argument-principle certificate that a small box
B = [x-r, x+r] x [y-r, y+r] contains exactly k zeros of Lambda_Q(s) (hence of Z_Q(s), since the
Gamma factor has no zeros and the poles 0,1 lie outside B).
Method: split dB into N segments; on each segment S evaluate Lambda_Q on a complex ball enclosing S
(theta.py gives a rigorous enclosure incl. tail bound).  If the resulting rectangle W excludes 0 and
lies in an open half-plane through 0 (checked as |Re mid|>rad_re or |Im mid|>rad_im), the continuous
arg change along S equals Arg(F(end)/F(start)) computed in balls. The total /2pi must enclose one integer.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json
sys.path.insert(0,'.')
from theta import Lambda
from flint import acb, arb, ctx

def ball_seg(Q,p,q,bits):
    m=(p+q)/2; hr=abs((q-p).real)/2; hi=abs((q-p).imag)/2
    s=acb(arb(m.real,hr),arb(m.imag,hi))
    return Lambda(Q,s,bits=bits)

def half_plane_ok(W):
    re,im=W.real,W.imag
    return (re>0) or (re<0) or (im>0) or (im<0)

def certify(Q,x,y,r,N=16,bits=80,depth=0):
    ctx.prec=bits+40
    c=[complex(x-r,y-r),complex(x+r,y-r),complex(x+r,y+r),complex(x-r,y+r)]
    pts=[]
    for i in range(4):
        a,b=c[i],c[(i+1)%4]
        pts+= [a+(b-a)*k/N for k in range(N)]
    pts.append(pts[0])
    tot=arb(0); worst=0
    for i in range(len(pts)-1):
        p,q=pts[i],pts[i+1]
        W=ball_seg(Q,p,q,bits)
        if not half_plane_ok(W):
            return None
        fp=Lambda(Q,acb(p.real,p.imag),bits=bits); fq=Lambda(Q,acb(q.real,q.imag),bits=bits)
        tot+= (fq/fp).arg()   # principal arg; valid since image of S lies in a half-plane
    k=tot/(2*arb.pi())
    lo=int(float(k.mid())+0.5)
    if (k-lo).abs_upper()<0.5:   # encloses exactly one integer
        return lo, str(k)
    return None

def run(Q,x,y,r=1e-6):
    for N in [8,16,32,64,128]:
        out=certify(Q,x,y,r,N=N)
        if out is not None:
            return dict(Q=Q,center=[x,y],halfwidth=r,segments_per_edge=N,count=out[0],count_ball=out[1])
    return dict(Q=Q,center=[x,y],halfwidth=r,count=None)

if __name__=='__main__':
    Q=tuple(int(a) for a in sys.argv[1:4]); x=float(sys.argv[4]); y=float(sys.argv[5])
    print(run(Q,x,y))
