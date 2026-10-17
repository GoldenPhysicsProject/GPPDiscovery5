"""Grok S57: independent check of the S14 real-space archimedean formula against the Fourier/digamma form
arch(g) = (1/2pi) int_R ghat(t) [Re psi(1/4+it/2) - log pi] dt  for Gaussians g(x)=exp(-x^2/(2s^2)). mpmath only."""
import mpmath as mp
mp.mp.dps = 40
def real_arch(g, g0):
    k1 = lambda x: (g0 - g(x)) * mp.exp(-x) / mp.sinh(x)
    k2 = lambda x: g(x) * mp.exp(-x / 4) * mp.sinh(3 * x / 4) / mp.sinh(x)
    return (-mp.euler - mp.log(mp.pi)) * g0 + mp.quad(k1, [0, 1, 5, 20, mp.inf]) - 2 * mp.quad(k2, [0, 1, 5, 20, mp.inf])
for s in [mp.mpf('0.3'), mp.mpf(1), mp.mpf(2)]:
    g = lambda x: mp.exp(-x**2 / (2 * s**2)); gh = lambda t: s * mp.sqrt(2 * mp.pi) * mp.exp(-s**2 * t**2 / 2)
    four = mp.quad(lambda t: gh(t) * (mp.re(mp.digamma(mp.mpf(1)/4 + 1j * t / 2)) - mp.log(mp.pi)), [0, 1, 5, 20, mp.inf]) / mp.pi
    r = real_arch(g, 1)
    print('s=%s real=%s fourier=%s diff=%s' % (s, mp.nstr(r, 25), mp.nstr(four, 25), mp.nstr(r - four, 5)))
# full explicit-formula normalization check: g=exp(-x^2/2), sum over zeros 2*sqrt(2pi)*sum exp(-gamma^2/2) ~ 1e-43
mp.mp.dps = 30
def primerange(a,b):
    sv=bytearray([1])*b; sv[0]=sv[1]=0
    for i in range(2,int(b**.5)+1):
        if sv[i]: sv[i*i::i]=bytearray(len(sv[i*i::i]))
    return [i for i in range(a,b) if sv[i]]
g = lambda x: mp.exp(-x**2 / 2)
pole = 2 * mp.sqrt(2 * mp.pi) * mp.exp(mp.mpf(1) / 8)
pr = mp.mpf(0)
for q in primerange(2, 600000):
    pk = q
    while pk < 600000:
        pr += 2 * mp.log(q) / mp.sqrt(pk) * g(mp.log(pk)); pk *= q
zs = 2 * mp.sqrt(2 * mp.pi) * sum(mp.exp(-mp.zetazero(k).imag**2 / 2) for k in range(1, 4))
print('pole - primes + arch =', mp.nstr(pole - pr + real_arch(g, 1), 15), ' zero side =', mp.nstr(zs, 5))
