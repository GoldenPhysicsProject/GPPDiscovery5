"""C.17: Thm C.15.A-type Schur certificate at the right endpoint 2L = log 3 (general-L / prime-set code c16_schur_general),
both parities; loss decomposition vs full-Q Galerkin (Ritz upper bounds). Float, NUMERICAL.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
from c16_schur_general import *
L=np.log(3)/2; P=prime_set(L); print("2L = log 3, P =",[(n,float(c)) for n,c in P],"(n = 3 at the endpoint contributes exactly 0)")
G={}
for par in ('even','odd'):
    for N in (128,256): G[(par,N)]=galerkin_full(L,P,N,par)
    print(f"{par}: full-Q Galerkin (Ritz upper bound) N=128 {G[(par,128)]:.4e}, N=256 {G[(par,256)]:.4e}",flush=True)
for par in ('even','odd'):
    for T,N in ((40,128),(60,128),(60,256),(100,256),(200,256)):
        for maj in ('analytic','sharp'):
            r=cert(L,P,float(T),N,par,maj); g=G[(par,N)]
            print(f"{par:4s} T={T:3d} N={N} {maj:8s}: lam_cert={r['lam']:+.3e} | Galerkin {g:.4e} = mu_A {r['muA']:+.3e} + comparison loss {g-r['muA']:.2e}; "
                  f"Schur loss mu_A-lam {r['muA']-r['lam']:.2e} (beta^2={r['B']**2:.2e}, c_tail={r['ctail']:.4f}, s_V={r['sV']:.3f}, eps^abs={r['eA']:.2e})",flush=True)
