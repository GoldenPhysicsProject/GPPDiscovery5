"""Round 11c: O7 sharpened. Low-frequency test of the explicit formula with Gaussian test functions
hhat(t) = exp(-t^2/(2 s^2)) (spectral mass concentrated in |t| < ~3s, i.e. inside the classical zero-free band |t| < 14.13).
Terms in the in-tree normalisation (same as K = pole - prime + arch):
  pole  = hhat(i/2)+hhat(-i/2) = 2 exp(1/(8 s^2))
  prime = sum_n Lambda(n) n^{-1/2} (h(log n)+h(-log n)),  h(x) = (s/sqrt(2pi)) exp(-s^2 x^2/2)
  arch  = (1/2pi) int hhat(t) [Re digamma(1/4+it/2) - log pi] dt   (split at t0 = 6.2898: arch_- , arch_+)
  W = pole - prime + arch  vs  zero side Z = sum_gamma hhat(gamma) (comparison only; zeros read from Daniel's
  Platt table zeros_14.dat, cross-checked against mpmath.zetazero).  NUMERICAL test of a step; no zeros as input."""
import mpmath as mp, struct, json, math
mp.mp.dps = 40
# --- read Platt block format (Bober reader logic) for comparison zeros
def platt(fn, maxz=3000):
    out = []
    with open(fn, 'rb') as f:
        nb = struct.unpack('Q', f.read(8))[0]
        for b in range(nb):
            hdr = f.read(32)
            if len(hdr) < 32: break
            t0, t1, N0, N1 = struct.unpack('ddQQ', hdr)
            Z = 0; t0m = mp.mpf(t0)
            for k in range(N1 - N0):
                z1, z2, z3 = struct.unpack('QIB', f.read(13))
                Z += (z3 << 96) + (z2 << 64) + z1
                out.append(t0m + mp.mpf(Z) * mp.mpf(2)**-101)
                if len(out) >= maxz: return out
    return out
gam = platt('/workspace/zeros/zeros_14.dat')
chk = [float(abs(gam[k] - mp.zetazero(k+1).imag)) for k in (0, 1, 9, 99)]
# --- primes
NMAX = 200000
sieve = bytearray([1])*(NMAX+1); sieve[0:2] = b'\x00\x00'
for i in range(2, int(NMAX**0.5)+1):
    if sieve[i]: sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
pps = []
for p in range(2, NMAX+1):
    if sieve[p]:
        q = p
        while q <= NMAX: pps.append((q, p)); q *= p
w = lambda t: mp.re(mp.digamma(mp.mpf(1)/4 + 1j*t/2)) - mp.log(mp.pi)
t0 = mp.findroot(w, 6.29)
rows = []
for s in ('1', '1.5', '2', '2.5', '3', '4'):
    s = mp.mpf(s)
    hh = lambda t: mp.exp(-t*t/(2*s*s))
    h = lambda x: s/mp.sqrt(2*mp.pi)*mp.exp(-s*s*x*x/2)
    pole = 2*mp.exp(1/(8*s*s))
    prime = mp.fsum(mp.log(p)/mp.sqrt(n)*2*h(mp.log(n)) for n, p in pps if s*mp.log(n) < 14)
    am = mp.quad(lambda t: hh(t)*w(t), [0, t0])/mp.pi
    ap = mp.quad(lambda t: hh(t)*w(t), [t0, 20*s, 60*s])/mp.pi
    W = pole - prime + am + ap
    Z = 2*mp.fsum(hh(g) for g in gam)
    rows.append(dict(s=float(s), pole=mp.nstr(pole, 12), prime=mp.nstr(prime, 12), arch_minus=mp.nstr(am, 12), arch_plus=mp.nstr(ap, 12),
                     W=mp.nstr(W, 6), zero_side=mp.nstr(Z, 6), W_minus_Z=mp.nstr(W - Z, 3),
                     lowband_bracket_arch_minus_prime_pole=mp.nstr(am - prime + pole, 12), max_term=mp.nstr(max(abs(pole), abs(prime), abs(am), abs(ap)), 4)))
    print(rows[-1], flush=True)
json.dump(dict(label='NUMERICAL', t0=mp.nstr(t0, 10), platt_vs_mpmath_absdiff=chk, n_zeros=len(gam), gamma_max=mp.nstr(gam[-1], 8), rows=rows),
          open('o7_lowband_results.json', 'w'), indent=1)
print('t0', t0, 'platt check', chk, len(gam), gam[-1])
