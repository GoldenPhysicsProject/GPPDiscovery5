"""Round 17, part B: Taylor coefficients of Xi_M from moments of the BPY law (integer data only), Turan and
Jensen-polynomial hyperbolicity (GORZ normalisation gamma(n) = n! m_{2n}/(2n)!, up to a positive constant),
  J^{d,n}(X) = sum_j C(d,j) gamma(n+j) X^j,   RH_M  <=>  all J^{d,n} hyperbolic (Polya 1927).
Root counting: python-flint fmpz_poly.complex_roots on the exact binary values of the 230-digit coefficients
(rigorous isolation for the rounded polynomial; coefficient error ~10^-200 relative each). Output r17_jensen.json"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from flint import fmpz_poly, fmpz
from laws17 import *
mp.mp.dps = 230
P = 200
H = mp.mpf('0.01')
K = 320
N_SHIFT = 150
D_HIGH = [10, 20, 40, 60, 80, 100, 140, 180, 240, 300]

def moments(name):
    if name in ('ZETA', 'EPS', 'DH'):
        us, ph = grid(Model(name), H)
        return moments_even(us, ph, H, K)
    us, ph = grid(Model('ZETA'), H)
    if name == 'FTH':
        th = mp.mpf('0.2'); Ma = tilted_moments(us, ph, H, 2*K, -th); Mb = tilted_moments(us, ph, H, 2*K, th)
    else:
        tau = mp.mpf(3); Ma = tilted_moments(us, ph, H, 2*K, -1j*tau); Mb = tilted_moments(us, ph, H, 2*K, 1j*tau)
    return [mp.re(x) for x in product_moments(Ma, Mb, K)]

def nreal(coeffs):
    # exact binary conversion of every coefficient (no common-scale truncation of small coefficients)
    parts = [mp.mpf(c).man_exp for c in coeffs]
    e0 = min(e for _, e in parts)
    ints = [fmpz(int(m_) << (e - e0)) for m_, e in parts]
    while ints and ints[-1] == 0: ints.pop()
    p = fmpz_poly(ints); deg = p.degree()
    roots = p.complex_roots()
    nr = sum(mult for r, mult in roots if r.imag == 0)
    return nr, deg, roots

def jensen(gam, d, n):
    return [mp.binomial(d, j)*gam[n + j] for j in range(d + 1)]

out = {}
for name in ['ZETA', 'EPS', 'DH', 'FTH', 'FTAU']:
    t0 = time.time()
    m = moments(name)
    gam = [mp.factorial(k)*m[k]/mp.factorial(2*k) for k in range(K + 1)]
    rec = dict(m0=mp.nstr(m[0], 12), all_moments_positive=all(x > 0 for x in m))
    # Turan: gamma(n)^2 - gamma(n-1) gamma(n+1) >= 0  (d=2)
    tur = [(gam[n]**2 - gam[n-1]*gam[n+1])/gam[n]**2 for n in range(1, K)]
    bad = [n + 1 for n, v in enumerate(tur) if v < 0]
    rec['turan_d2'] = dict(n_range=[1, K - 1], failures=bad[:20], n_fail=len(bad), min_rel=mp.nstr(min(tur), 6))
    # low degree, all shifts
    low = {}
    for d in range(3, 9):
        fails = []
        for n in range(0, N_SHIFT + 1):
            nr, deg, _ = nreal(jensen(gam, d, n))
            if nr < deg: fails.append(n)
        low[str(d)] = dict(n_range=[0, N_SHIFT], non_hyperbolic_n=fails[:30], n_fail=len(fails))
    rec['jensen_low_degree'] = low
    # high degree, small shifts
    high = {}
    for n in [0, 1, 2, 5, 10]:
        for d in D_HIGH:
            if n + d > K: continue
            nr, deg, roots = nreal(jensen(gam, d, n))
            high[f'd{d}_n{n}'] = dict(d=d, n=n, real_roots=int(nr), degree=int(deg), hyperbolic=bool(nr == deg))
    rec['jensen_high_degree'] = high
    first = {}
    for n in [0, 1, 2, 5, 10]:
        ds = [v['d'] for k, v in high.items() if v['n'] == n and not v['hyperbolic']]
        first[str(n)] = min(ds) if ds else None
    rec['first_non_hyperbolic_d_by_shift'] = first
    out[name] = rec
    print(name, round(time.time() - t0, 1), 's', rec['all_moments_positive'], rec['turan_d2']['n_fail'],
          {d: v['n_fail'] for d, v in low.items()}, first, flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r17_jensen.json'), 'w'), indent=1)
