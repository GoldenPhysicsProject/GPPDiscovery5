"""C.16 odd at the p=2 right endpoint 2L = log 3 (ell = log 1.5): Thm C.16.A float lambda_cert (NUMERICAL), T/N scan.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
from c16_schur_general import *
L=np.log(3)/2; P=prime_set(L); print("2L = log 3, P =",[(n,float(c)) for n,c in P])
for N in (64,128,256): print(f"full-Q odd Galerkin N={N} (upper bound): {galerkin_full(L,P,N,'odd'):+.6e}",flush=True)
for T,N in ((40,96),(60,64),(60,128),(60,192),(80,128),(80,192),(100,192),(100,256)):
    for maj in ('analytic','sharp'):
        r=cert(L,P,float(T),N,'odd',maj)
        print(f"{maj:8s}",fmt(r),f"| comparison loss (Gal N - mu_A) {galerkin_full(L,P,N,'odd')-r['muA']:.3e} | Schur loss {r['muA']-r['lam']:.3e}",flush=True)
