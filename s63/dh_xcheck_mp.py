"""Grok S63: high-precision (mpmath) cross-check of a DH Ritz witness, using the Claude/Codex DH explicit formula
(odd-gamma kernel e^{-3y/2}/(1-e^{-2y}), constant gamma+log(pi/5)+log(1-e^{-4L}), Dirichlet-inverse c(n)) instead of
dhritz.py's S47 sech-difference arch term. Autocorrelation via the closed form (selftest.py, 4e-31), reduced to O(K) per point:
  F(t) = sum_m v_m^2 F_mm(t) + 2 sum_n C_n sin(a_n t),  C_n = v_n sum_{m!=n} v_m eps_mn a_m/(a_m^2-a_n^2).
NUMERICAL (tanh-sinh on 2*npan panels), not Arb. Usage: dh_xcheck_mp.py cosvec.json [dps] [npan]"""
import sys, json, mpmath as mp
from fractions import Fraction
d=json.load(open(sys.argv[1])); mp.mp.dps=int(sys.argv[2]) if len(sys.argv)>2 else 50; npan=int(sys.argv[3]) if len(sys.argv)>3 else 120
L=mp.mpf(Fraction(d["L"]).numerator)/Fraction(d["L"]).denominator; K=d["K"]; sec=d["sector"]
v=[mp.mpf(Fraction(s).numerator)/Fraction(s).denominator for s in d["vec"]]
a=[(2*m+1)*mp.pi/(2*L) for m in range(K)] if sec=="even" else [m*mp.pi/L for m in range(1,K+1)]
C=[v[n]*mp.fsum(v[m]*(-1)**(m+n)*a[m]/(a[m]**2-a[n]**2) for m in range(K) if m!=n) for n in range(K)]
def F(t): return mp.fsum(v[m]**2*((2*L-t)*mp.cos(a[m]*t)/2+mp.sin(a[m]*t)/(2*a[m])) + 2*C[m]*mp.sin(a[m]*t) for m in range(K))
F0=F(mp.mpf(0)); s5=mp.sqrt(5); kap=(mp.sqrt(10-2*s5)-2)/(s5-1); ach={0:0,1:1,2:kap,3:-kap,4:-1}
N=int(mp.floor(mp.e**(2*L))); c=[mp.mpf(0)]*(N+1)
for n in range(2,N+1): c[n]=ach[n%5]*mp.log(n)-mp.fsum(ach[q%5]*c[n//q] for q in range(2,n) if n%q==0)
arith=2*mp.fsum(c[n]/mp.sqrt(n)*F(mp.log(n)) for n in range(2,N+1) if mp.log(n)<2*L)
def kern(y): return (mp.e**(-1.5*y)*F(y)-mp.e**(-2*y)*F0)/(-mp.expm1(-2*y))
I=mp.quad(kern, mp.linspace(0,2*L,npan+1))
arch=(mp.euler+mp.log(mp.pi/5)+mp.log(-mp.expm1(-4*L)))*F0+2*I
Q=-(arch+arith)
print(json.dumps({"file":sys.argv[1],"sector":sec,"L":d["L"],"K":K,"dps":mp.mp.dps,"npan":npan,"F0":mp.nstr(F0,20),"LV2":mp.nstr(L*mp.fsum(x*x for x in v),20),
 "Q_over_norm":mp.nstr(Q/F0,15),"dhritz_arb":d["RQ"]}))
