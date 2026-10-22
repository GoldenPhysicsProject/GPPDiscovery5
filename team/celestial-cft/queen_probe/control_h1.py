"""Control for epstein_sinh.py: the same Epstein-normalised window form for x^2+y^2 (h=1, Z = 4 zeta(s) L(s,chi_-4), Euler product, GRH expected)
on the same vectors. If the normalisation (Omega_E, nu_E, pole) is right, Q should NOT run away to -infinity. NUMERICAL.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys; sys.argv=['x']
import importlib.util, mpmath as mp, json
spec=importlib.util.spec_from_file_location('es','epstein_sinh_lib.py'); es=importlib.util.module_from_spec(spec); spec.loader.exec_module(es)
cE=es.logderiv(es.rQ((1,0,1))/4.0*0+ (lambda r:(r/r[1]))(es.rQ((1,0,1))))
cE[0]=0
Om=float(mp.log(4/(4*mp.pi**2))+2*mp.digamma(0.5))
out=[]
for (x0,y0) in [(0.9079933666439008,89.55982797855877),(0.9329696974854141,15.668249531278473)]:
    for L in [2.0,3.0,4.0,5.0,6.0]:
        r=es.run(L,y0,x0-0.5,cE,Om,lambda u:1/mp.sinh(u/2),True); out.append(dict(rho0=[x0,y0],**r))
        print(f"x^2+y^2 on vector rho0={x0:.4f}+{y0:.3f}i L={L}: Q/|f|^2={r['Q_over_norm']:+.4g} (arch {r['arch']:+.3f}, prime {r['prime']:+.3f})",flush=True)
json.dump(out,open('control_h1.json','w'),indent=1)
