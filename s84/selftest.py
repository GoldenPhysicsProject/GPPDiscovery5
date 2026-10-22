# check the fhat closed forms by quadrature, and the explicit formula normalisation on a smooth-ish low-L vector
import mpmath as mp, json
from spec_decomp import fhat_factory
mp.mp.dps=30
for sec in ['even','odd']:
    vec=['1','-1/3','1/7','2/5']; f,n2,L,v,fr,_=fhat_factory(vec,'7/10',sec)
    for t in [3.3, 14.1347, 40.0]:
        g=(lambda u: sum(v[m]*mp.cos(fr[m]*u) for m in range(4))) if sec=='even' else (lambda u: sum(v[m]*mp.sin(fr[m]*u) for m in range(4)))
        q=mp.quad(lambda u:g(u)*mp.exp(1j*t*u),[-L,0,L])
        print(sec,t,mp.nstr(abs(q),15),mp.nstr(abs(f(t)),15))
