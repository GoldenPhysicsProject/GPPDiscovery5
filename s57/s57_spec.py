"""Grok S57 H57c: spectral anatomy of the Weil minimizer (odd vs even). Legendre basis (EDGEK=0), f(x)=sum c_j P_j(x/L):
fhat(t) = 2L sum_j c_j i^j j_j(Lt) (spherical Bessel, exact). Frequency form of the Weil functional:
  Q(f) = pole(f) + (1/2pi) int |fhat(t)|^2 [Omega(t) - Pi_L(t)] dt,  Omega = Re psi(1/4+it/2) - log pi,
  Pi_L(t) = sum_{log n<2L} (2 Lambda(n)/sqrt n) cos(t log n);  pole = 2 fhat(i/2)fhat(-i/2) (<0 for odd f).
For the Ritz minimizer we report: band masses, where the arch is negative (Omega<0, |t|<t_Omega), where the primes lift
(Pi_L<0), spectral centroid, and the reconstructed Q (consistency vs Ritz lambda). NUMERICAL. RH OPEN.
usage: EDGEK=0 python s57_spec.py SECTOR M L1,L2,..."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's56')); sys.path.insert(1, os.path.join(HERE, '..', 's53'))
os.chdir(os.path.join(HERE, '..', 's56'))
assert os.environ.get('EDGEK') == '0'
_a = sys.argv; sys.argv = sys.argv[:1]
from s53_odd_floor import build_all, tomp, weight, vonmangoldt
sys.argv = _a
import mpmath as mp
mp.mp.dps = 60
sec, M = sys.argv[1], int(sys.argv[2])
Om = lambda t: mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * t / 2)) - mp.log(mp.pi)
tO = mp.findroot(Om, 6.0)
for Ls in sys.argv[3].split(','):
    t0 = time.time(); L = mp.mpf(Ls); I = list(range(0 if sec == 'even' else 1, M, 2))
    mp.mp.dps = 130
    A, P, T, ns = build_all(M, Ls, 1024)
    W = tomp(A, I) + tomp(P, I)
    for n in ns: W -= weight(n) * tomp(T[n], I)
    Gd = [2 * L / (2 * j + 1) for j in I]   # Legendre Gram is diagonal
    S = mp.matrix(len(I), len(I))
    for a in range(len(I)):
        for b in range(len(I)): S[a, b] = W[a, b] / mp.sqrt(Gd[a] * Gd[b])
    ev, Q = mp.eigsy(S); i0 = min(range(len(I)), key=lambda i: ev[i])
    c = [Q[a, i0] / mp.sqrt(Gd[a]) for a in range(len(I))]   # ||f||=1
    import numpy as np
    from scipy.special import spherical_jn, digamma
    cf = np.array([float(x) for x in c]); Lf = float(L); tmax = 60 + 4 * M / Lf
    t = np.linspace(0, tmax, 400001); dt = t[1] - t[0]
    fh = np.zeros_like(t, dtype=complex)
    for a, j in enumerate(I): fh += 2 * Lf * cf[a] * (1j)**j * spherical_jn(j, Lf * t)
    dens = np.abs(fh)**2 / (2 * np.pi)
    Omv = np.real(digamma(0.25 + 0.5j * t)) - np.log(np.pi)
    Piv = sum(float(weight(n)) * np.cos(t * np.log(n)) for n in ns)
    w = np.full_like(t, dt); w[0] = w[-1] = dt / 2; w *= 2   # even integrand: int_R = 2 int_0^inf
    mass = (w * dens).sum(); archQ = (w * dens * Omv).sum(); primQ = -(w * dens * Piv).sum()
    # pole: fhat(i/2) fhat(-i/2) with f odd/even real: compute directly in mpmath
    def fhz(z):
        return 2 * L * mp.fsum(c[a] * (1j)**j * mp.sqrt(mp.pi / (2 * L * z)) * mp.besselj(j + 0.5, L * z) for a, j in enumerate(I))
    pole = 2 * float(mp.re(fhz(1j / 2) * fhz(-1j / 2)))   # ghat(i/2)+ghat(-i/2)
    neg = t < float(tO); lift = np.maximum(-Piv, 0); lob = (t > np.pi / (2 * np.log(2))) & (t < 3 * np.pi / (2 * np.log(2)))
    r = dict(sec=sec, M=M, L=Ls, ns=ns, lam_ritz=mp.nstr(ev[i0], 8), Q_rec=float(pole + archQ + primQ), mass=float(mass), pole=pole,
        arch=float(archQ), primes=float(primQ), t_Omega=float(tO), mass_frac_Omega_neg=float((w * dens * neg).sum() / mass),
        mass_frac_lobe2=float((w * dens * lob).sum() / mass), mass_frac_t_lt_1=float((w * dens * (t < 1)).sum() / mass),
        arch_from_Omega_neg=float((w * dens * Omv * neg).sum()), arch_from_Omega_pos=float((w * dens * Omv * ~neg).sum()),
        prime_lift_pos=float((w * dens * lift).sum()), prime_drag=float(-(w * dens * np.maximum(Piv, 0)).sum()),
        lift_frac_in_Omega_neg=float((w * dens * lift * neg).sum() / (w * dens * lift).sum()),
        centroid=float((w * dens * t).sum() / mass), t_peak=float(t[np.argmax(dens)]), wall_s=round(time.time() - t0))
    g1 = 14.134725141734695; r['mass_frac_below_gamma1'] = float((w * dens * (t < g1)).sum() / mass)
    if os.environ.get('ZS'):
        zs = json.load(open(os.path.join(HERE, 'spec', 'zeros600.json'))); mp.mp.dps = 40; acc = mp.mpf(0); part = {}
        for k, gz in enumerate(zs):
            acc += 2 * abs(fhz(mp.mpf(gz)))**2   # Q = sum over zeros (+-gamma) of ghat(gamma), ghat=|fhat|^2
            if k + 1 in (1, 5, 20, 100, 600): part[k + 1] = mp.nstr(acc, 8)
        r['zero_side_partial'] = part   # partial zero sums, first N zeros (RH verified there)
    print(json.dumps(r), flush=True)
# Grok S57 (2026-10-17). RH OPEN. phi absent unless reported.
