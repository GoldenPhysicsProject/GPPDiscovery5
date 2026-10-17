"""Grok S63: INDEPENDENT float64 cross-check of a DH Ritz witness from dhritz.py (FORM=DH, S47 sech-difference arch).
Evaluates the Claude/Codex DH explicit formula (work/claude-code/2026-10-13_dh_modulated_window, CODEX 2026-10-13 witness note):
  Q_DH(g) = -(gamma+log(pi/5)+log(1-e^{-4L}))G0 - 2 int_0^{2L} [e^{-3y/2}G(y)-e^{-2y}G0]/(1-e^{-2y}) dy - 2 sum_{n<e^{2L}} c(n) n^{-1/2} G(log n)
with c(n) the -D'/D Dirichlet coefficients (recurrence c(n)=a(n)log n - sum_{d|n,1<d<n} a(d)c(n/d)) recomputed here in floats,
and G(x)=int f(y)f(y-x)dy computed by DIRECT Gauss-Legendre quadrature of the trig expansion (no closed-form autocorrelation).
Usage: dh_xcheck.py cosvec_DH_<sec>_L.._K...json"""
import sys, json, numpy as np
from fractions import Fraction
d=json.load(open(sys.argv[1])); L=float(Fraction(d["L"])); K=d["K"]; sec=d["sector"]
v=np.array([float(Fraction(s)) for s in d["vec"]])
if sec=="even": fr=np.array([(2*m+1)*np.pi/(2*L) for m in range(K)]); bas=np.cos
else: fr=np.array([m*np.pi/L for m in range(1,K+1)]); bas=np.sin
def f(y): return bas(np.outer(y,fr))@v
xg,wg=np.polynomial.legendre.leggauss(1600)
def G(x):  # int_{x-L}^{L} f(y) f(y-x) dy
    a,b=x-L,L; y=(b-a)/2*xg+(a+b)/2; return (b-a)/2*np.sum(wg*f(y)*f(y-x))
G0=G(0.0)
s5=5**.5; kap=((10-2*s5)**.5-2)/(s5-1)
a={0:0,1:1,2:kap,3:-kap,4:-1}
N=int(np.exp(2*L))+1; c=[0.0]*(N+1)
for n in range(2,N+1):
    c[n]=a[n%5]*np.log(n)-sum(a[dd%5]*c[n//dd] for dd in range(2,n) if n%dd==0)
arith=2*sum(c[n]/np.sqrt(n)*G(np.log(n)) for n in range(2,N+1) if np.log(n)<2*L)
# arch integral, split in 64 panels x 64 GL nodes
xp,wp=np.polynomial.legendre.leggauss(64); edges=np.linspace(0,2*L,65); I=0.0
for lo,hi in zip(edges[:-1],edges[1:]):
    for t,w in zip((hi-lo)/2*xp+(lo+hi)/2,(hi-lo)/2*wp):
        I+=w*(np.exp(-1.5*t)*G(t)-np.exp(-2*t)*G0)/(-np.expm1(-2*t))
arch=(np.euler_gamma+np.log(np.pi/5)+np.log(-np.expm1(-4*L)))*G0+2*I
Q=-(arch+arith)
print(json.dumps({"file":sys.argv[1],"L":d["L"],"sector":sec,"K":K,"G0":G0,"L*|v|^2":L*float(v@v),"Q_over_norm":Q/G0,"arb_RQ_dhritz":d["RQ"]}))
