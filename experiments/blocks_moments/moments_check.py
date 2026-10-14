"""High-precision verification suite for modular thermality of P(lambda)."""
import mpmath as mp
mp.mp.dps=60
P=lambda l: mp.pi*l/mp.sinh(mp.pi*l) if l else mp.mpf(1)
phat=lambda x: mp.mpf(1)/(4*mp.cosh(x/2)**2)
def rel(a,b): return abs(a-b)/max(1,abs(a),abs(b))
# Fourier pair with the convention F(lambda)=int e^{-i lambda x} f(x) dx.
for l in [mp.mpf('0'),mp.mpf('.3'),mp.mpf('1.1'),mp.mpf('2.7')]:
    num=mp.quad(lambda x:phat(x)*mp.cos(l*x),[-mp.inf,0,mp.inf])
    assert rel(num,P(l))<mp.mpf('1e-40')
# Planck representation P=2*pi*lambda*(n_B(pi lambda)-n_B(2 pi lambda)).
for l in [mp.mpf('.1'),mp.mpf('.7'),mp.mpf('2.3')]:
    nb=lambda q:1/mp.expm1(q)
    rhs=2*mp.pi*l*(nb(mp.pi*l)-nb(2*mp.pi*l))
    assert rel(rhs,P(l))<mp.mpf('1e-45')
# Weierstrass/log expansion derivative check: d/dlambda log P = 1/lambda - pi*coth(pi lambda).
for l in [mp.mpf('.2'),mp.mpf('1.1'),mp.mpf('3')]:
    lhs=mp.diff(P,l)/P(l)
    rhs=1/l-mp.pi*mp.coth(mp.pi*l)
    assert rel(lhs,rhs)<mp.mpf('1e-45')
# Matsubara poles: residue of P at lambda=i n is i*(-1)^n*n for n != 0.
for n in [1,2,3]:
    z=1j*n
    # residue via (lambda-z)P(lambda), evaluated close to z.
    eps=mp.mpf('1e-20')
    res=(eps)*P(z+eps)
    exact=1j*((-1)**n)*n
    assert rel(res,exact)<mp.mpf('1e-15')
print('PASS: Fourier pair P <-> sech^2(x/2)/4')
print('PASS: Planck representation at beta=2*pi')
print('PASS: logarithmic derivative / Weierstrass-compatible identity')
print('PASS: Matsubara pole residues')
