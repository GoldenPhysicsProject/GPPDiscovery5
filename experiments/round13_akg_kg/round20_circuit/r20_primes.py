"""Task 3: arithmetic inside the ladder. Re-centre at u0 = 100 (s0 = 1/2 + 10 = 10.5, inside Re s > 1) where the
explicit formula holds with absolutely convergent prime sums:
  m(u) = [1/s + 1/(s-1) - (1/2) log pi + (1/2) psi(s/2) - sum_n Lambda(n) n^{-s}] / (2 sqrt u),  s = 1/2 + sqrt u.
Moments d_j: m(u0+v) = sum_j (-1)^j d_j v^j (= sum_rho (t_rho+u0)^{-j-1} under RH). Prime part of d_j:
  -sum_n Lambda(n) n^{-1/2} (1/j!) (-d/du)^j [e^{-sqrt(u) log n}/(2 sqrt u)]  =  -sum_n Lambda(n) n^{-1/2} int_0^inf tau^j/j! e^{-u0 tau}
  e^{-(log n)^2/(4 tau)} / (2 sqrt(pi tau)) dtau  < 0  for every j and every prime power (PROVED, Laplace form).
Route A: Taylor coefficients of xi at 1/2 (r20_core), shifted to u0. Route B: arch + primes <= P. Then S-fraction at u0."""
import json, math
import mpmath as mp
from r20_core import F_coeffs, series_inv, m_series, s_fraction, ladder
U0 = 100; J = 40
# ---- Route A
K = 260; D = 500
f = F_coeffs('ZETA', K, D); mp.mp.dps = D
u0 = mp.mpf(U0)
Fs = [mp.fsum(f[k]*mp.binomial(k, i)*u0**(k - i) for k in range(i, K + 1)) for i in range(J + 2)]
Fs2 = [mp.fsum(f[k]*mp.binomial(k, i)*u0**(k - i) for k in range(i, K - 39)) for i in range(J + 2)]
conv = float(max(abs(a - b)/abs(a) for a, b in zip(Fs, Fs2)))
muA = m_series(Fs, J)
# ---- Route B
mp.mp.dps = 80
u0 = mp.mpf(U0); r0 = mp.sqrt(u0)
def arch(v):
    u = u0 + v; s = mp.mpf(1)/2 + mp.sqrt(u)
    return (1/s + 1/(s - 1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2)/(2*mp.sqrt(u))
archT = mp.taylor(arch, 0, J - 1)
# series of w = sqrt(u0+v) and 1/(2w)
w = [mp.binomial(mp.mpf(1)/2, i)*u0**(mp.mpf(1)/2 - i) for i in range(J)]
inv2w = [x/2 for x in series_inv(w, J)]
def exp_series(a, n):           # exp of series with a[0] = 0
    e = [mp.mpf(0)]*n; e[0] = mp.mpf(1)
    for k in range(1, n):
        e[k] = mp.fsum(j*a[j]*e[k - j] for j in range(1, k + 1))/k
    return e
def prime_series(P):
    tot = [mp.mpf(0)]*J
    ps = [p for p in range(2, P + 1) if all(p % q for q in range(2, int(p**0.5) + 1))]
    for p in ps:
        q = p
        while q <= P:
            Lq = mp.log(q); Lam = mp.log(p)
            e = exp_series([0] + [-Lq*w[i] for i in range(1, J)], J)
            pre = Lam*mp.e**(-Lq*(mp.mpf(1)/2 + r0))
            for k in range(J):
                tot[k] += pre*mp.fsum(e[i]*inv2w[k - i] for i in range(k + 1))
            q *= p
    return tot
out = dict(u0=U0, J=J, routeA_shift_convergence=conv)
mp.mp.dps = 80
d_A = [((-1)**j)*muA[j] for j in range(J)]
cA, aA = s_fraction(muA, J - 1)
out['d_A'] = [mp.nstr(x, 25) for x in d_A]
out['a_A'] = [mp.nstr(x, 20) for x in aA]
for P in (0, 2, 3, 5, 10, 30, 100, 1000, 10000):
    pr = prime_series(P) if P >= 2 else [mp.mpf(0)]*J
    muB = [a - b for a, b in zip(archT, pr)]
    d_B = [((-1)**j)*muB[j] for j in range(J)]
    rel = [float(abs(x - y)/abs(x)) for x, y in zip(d_A, d_B)]
    cB, aB = s_fraction(muB, J - 1)
    agree = 0
    while agree < min(len(aA), len(aB)) and abs(aB[agree] - aA[agree]) < 1e-6*abs(aA[agree]): agree += 1
    neg = next((k + 1 for k, x in enumerate(aB) if x <= 0), None)
    prime_part_signs = ''.join('-' if ((-1)**j)*(-pr[j]) < 0 else ('0' if pr[j] == 0 else '+') for j in range(min(J, 20)))
    row = dict(max_rel_err_moments_first10=max(rel[:10]), rel_err_moments=[float('%.3g' % x) for x in rel[:20]],
               n_ladder_elements_agreeing_1em6=agree, first_nonpositive_a=neg,
               a_first8=[mp.nstr(x, 8) for x in aB[:8]], sign_of_prime_contribution_to_d_j=prime_part_signs)
    out[f'P={P}'] = row
    print(P, row['max_rel_err_moments_first10'], agree, neg, prime_part_signs, flush=True)
json.dump(out, open('r20_primes.json', 'w'), indent=1)
