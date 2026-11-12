# Consistency check: full Weil form (all primes <= e^{2 sigma}) must be >= 0 if RH (known true at these heights in effect).
from wS_indefinite import run
import sympy
for sig in (0.3,1,2,3,4,5):
    S=list(sympy.primerange(2,int(__import__('math').exp(2*sig))+2))
    g0,pl,ai,W,lim=run(sig,S,n=16001)
    print(f"sigma={sig} #S={len(S)}: Ainf/g0={ai:.5f} W_full/g0={W:.6f}",flush=True)
