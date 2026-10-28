"""Grok audit S103i(2) of Celestial C.15 (NUMERICAL cross-checks; run from repo root). Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import numpy as np, sys
sys.path.insert(0,'team/celestial-cft/arp_c')
from scipy.special import digamma
from flint import arb, acb, ctx
ctx.prec=128
a=np.log(2)
def Om(t): return np.real(digamma(0.25+0.5j*np.asarray(t,dtype=complex)))-np.log(np.pi)
Om0=Om(0.0)
# K closed form check: colP[k] via direct Fourier integral (float), ell=a/6, M=28
M,ns=28,4; L=a*M/(2*(M-ns)); h=2*L/M
def colP_direct(k,X=4e5):
    x,w=np.polynomial.legendre.leggauss(40); tot=0.0
    edges=np.concatenate([np.linspace(0,200,401),np.geomspace(200,X,20000)])
    for i in range(len(edges)-1):
        t=0.5*(edges[i+1]-edges[i])*(x+1)+edges[i]; ww=0.5*(edges[i+1]-edges[i])*w
        f=(0.25+t*t)*(Om(t)-Om0)*h*h*np.sinc(t*h/2/np.pi)**4*np.cos(k*h*t)/np.pi
        tot+=np.sum(ww*f)
    avg={0:3/8,1:-1/4,2:1/16}.get(k,0.0)
    # tail: (1/pi) int_X^inf 16/h^2 * h^2/t^2*(Om(t)-Om0)*avg ; Om ~ log(t/2)-log pi
    tail=16*avg/(np.pi*h*h)*((np.log(X/2)-np.log(np.pi)-Om0)+1)/X
    return tot+tail
def K(d):
    d=abs(d)
    if d==0: return acb(2).zeta(acb(0.25)).real/4-acb(4).zeta(acb(0.25)).real/64
    s=arb(0)
    for n in range(4000):
        c=arb(2*n)+arb(0.5); s+=(-c*d).exp()*(1/c**2-1/(4*c**4))
    return s-d*acb(3).zeta(acb(0.25)).real/32
hA=2*(arb(2).log()*M/(2*(M-ns)))/M
Kv=[K(arb(k)*hA) for k in range(12)]
w=[1,-2,1]
for k in (0,1,2,3,5):
    s=sum(w[p]*w[q]*Kv[abs(k+p-q)] for p in range(3) for q in range(3))/hA**2
    print(k, float(s.mid()), colP_direct(k))
