"""Grok S61: spot check of exact claims in Codex PR #11 comments 6231548549 / 6231730580 (celestial W(1+inf), c=0 audit). NUMERICAL/exact."""
import mpmath as mp, json
from fractions import Fraction as Fr
mp.mp.dps = 40
out = {}
# (a) |Gamma(1+il)^3/Gamma(2+3il)|^2 closed form and limit
ds = []
for l in [mp.mpf('0.1'), mp.mpf('0.7'), mp.mpf(2), mp.mpf(5)]:
    lhs = abs(mp.gamma(1+1j*l)**3/mp.gamma(2+3j*l))**2
    rhs = mp.pi**2*l**2*(4+3*mp.csch(mp.pi*l)**2)/(3*(1+9*l**2))
    ds.append(float(abs(lhs-rhs)/rhs))
out['gamma_modsq_reldiff_max'] = max(ds)
out['gamma_modsq_limit_minus_4pi2_27'] = mp.nstr(mp.pi**2*mp.mpf(40)**2*(4+3*mp.csch(40*mp.pi)**2)/(3*(1+9*40**2)) - 4*mp.pi**2/27*(9*1600)/(1+9*1600), 3)
out['lambda0_limit'] = mp.nstr(mp.limit(lambda l: mp.pi**2*l**2*(4+3*mp.csch(mp.pi*l)**2)/(3*(1+9*l**2)), 0), 15)  # = 1/... (|Gamma(1)^3/Gamma(2)|^2 = 1)
# (b) Kac-Radul cocycle sum_{j=1}^n ((n+1)/2-j)^2 = n(n^2-1)/12, exact
out['kac_radul_ok_n_le_60'] = all(sum((Fr(n+1, 2)-j)**2 for j in range(1, n+1)) == Fr(n*(n*n-1), 12) for n in range(1, 61))
# (c) Virasoro [L_m,L_-m] = 2m L0 + c/12 (m^3-m): m=2 -> 4L0 + c/2; <0|L_2 L_-2|0> = c/2
out['vir_m2_central_coeff'] = str(Fr(2**3-2, 12))
# (d) shadow-cut weight w = l/(8 sinh(pi l)): int w = 1/16, Fourier = sech^2(t/2)/16; OS 2x2 witness
w = lambda l: l/(8*mp.sinh(mp.pi*l)) if l != 0 else 1/(8*mp.pi)
out['int_w_minus_1_16'] = mp.nstr(2*mp.quad(w, [0, mp.inf]) - mp.mpf(1)/16, 3)
for t in [mp.mpf('0.3'), mp.log(6)]:
    ft = 2*mp.quad(lambda l: w(l)*mp.cos(l*t), [0, mp.inf])
    out['fourier_minus_sech2_over16_t=%s' % mp.nstr(t, 5)] = mp.nstr(ft - mp.sech(t/2)**2/16, 3)
K = lambda e: Fr(4)*e/((e+1)**2)/16      # sech^2(t/2)/16 with e = e^t (rational at t = log of rational)
M = [[K(Fr(4)), K(Fr(6))], [K(Fr(6)), K(Fr(9))]]
out['OS_matrix'] = [[str(x) for x in r] for r in M]; out['OS_det'] = str(M[0][0]*M[1][1]-M[0][1]**2)
json.dump(out, open('celestial_spot.json', 'w'), indent=1); print(json.dumps(out, indent=1))
