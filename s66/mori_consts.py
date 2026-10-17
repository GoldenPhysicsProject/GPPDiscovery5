"""Grok S66: scalar constants of the Mori/Liu bounded comparison (Mori [9]=paper VII, zenodo 23265728:
Prop 3.1 tail, Lemma 3.2 projection errors e,n,h, alpha,beta of (3.3)), in Arb.  Validation vs Mori VII (log 12)
and Mori VIII (log 15).  Status: Arb ball evaluation of published formulas (CERTIFIED arithmetic, formulas Mori's)."""
import sys, json
from flint import arb, acb, ctx
ctx.prec = 256

def a_arch(t):
    return acb(arb(1)/4, arb(t)/2).digamma().real - arb.pi().log()

def consts(L, Om, s, m, u0, u1, N, rho=3, delta_log2=None):
    L = arb(L); Om = arb(Om); s = arb(s); m = arb(m); g = s + m; b = 2*s + m
    c = Om*L
    C = (Om/(2*arb.pi())).log() - 1/Om - g
    dc = (-2*(c+1)).exp()/(4*(c+1)*(c+4))
    E0 = (c+1)/(arb.pi()*c**2); E1 = (3*c**2+6*c+2)/(arb.pi()*c**3)
    rho = arb(rho)
    rN2 = 8*L*rho**2/(rho-1)**2 * (c*(rho-1/rho)).exp() * rho**(-2*N)
    rN = rN2.sqrt()
    K0 = Om/arb.pi() * max(abs(a_arch(0)-g).upper(), abs(a_arch(Om)-g).upper())
    Hp = (2*L).sqrt()*(L/2).cosh()
    qN = (2*L).sqrt()*(L/2).exp()*(L/2)**N/arb(N).fac()
    e = (K0*(2*L).sqrt() + (Om/arb.pi()).sqrt()*(arb(u0)*E0.sqrt()+arb(u1)*E1.sqrt()))*rN + 4*Hp*qN
    n = K0*rN2 + 4*qN**2
    h = (K0 + (arb(u0)+arb(u1))*Om/arb.pi())*rN2 + 4*qN**2
    al = e + 2*e**2/(m-n); be = (e+n)/m**2 + 2*h**2/(m**2*(m-n))
    out = dict(gamma=g, b=b, c=c, C=C, mlog2_dc=-dc.log()/arb(2).log(), u0E0=arb(u0)*E0, u1E1=arb(u1)*E1,
               e=e, n=n, h=h, alpha=al, beta=be, alpha_plus_beta_b2=al+be*b**2, K0=K0)
    ok = bool(C > 0) and bool(arb(u0)*E0 <= C) and bool(arb(u1)*E1 <= C) and bool(n < m) and bool(c > 2) and bool(Om >= arb(15)/4)
    return out, ok

def show(tag, out, ok):
    print(tag, 'hypotheses OK' if ok else 'HYPOTHESIS FAILS')
    for k, v in out.items(): print('  %-20s %s' % (k, v.str(8)))

if __name__ == '__main__':
    l12 = arb(12).log()/2; l15 = arb(15).log()/2
    show('Mori VII log12 (paper: C~0.12790, -log2 d=1514.27, e<2.08e-94, n<1.27e-191, h<3.82e-190, a+bb^2<2^-303)',
         *consts(l12, 416, arb(57)/16, arb(1)/2, 207, 68, 832))
    show('Mori VIII log15 (paper: C~0.136994, -log2 d=2900.29, u0E0~0.13684, u1E1~0.13634, e<2.75e-213, n<1.09e-429, h<6.29e-428, a+bb^2<2^-697)',
         *consts(l15, 736, arb(33)/8, arb(1)/2, 428, 142, 1664))
