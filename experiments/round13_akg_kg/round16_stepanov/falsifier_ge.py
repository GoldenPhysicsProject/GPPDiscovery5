"""Golden Eye AUDIT15 s6 falsifier for Round 16: run the Stepanov toy's OUTPUT unchanged on F_q-side fakes.
Setting = the Part 1 setting: probe-B elliptic curves E/F_p base-changed to F_q, q = p^2 (Bombieri needs q square > 16).
Data per model, n = 1..40 (all from the model's own L-numerator; exact integers where the model has integer data):
  CURVE  zeta_A of E/F_q        c_n = q^n - a_n,  a_n = alpha^n + conj(alpha)^n,  alpha = alpha_p^2
  Z0     probe-B class-0 partial zeta over F_q, Z_0 = (1 - qT + qT^2)/(1 - qT): FE, no Euler product, 1 ghost
  FTAU   unitary twist  F_{i tau}: c_n * 2cos(n tau log q)                (tau = 3, as on the Z side)
  FTH    non-unitary    F_theta : c_n * (q^{n theta} + q^{-n theta})      (theta = 0.2)
Errors e_n = c_n - main_n are computed exactly (integers: -a_n; twisted: -a_n * weight), no cancellation.
Toy inequalities (Round 16 Part 1, genus 1, constant 2g+1 = 3 per Frobenius piece; FTAU/FTH have 2 pieces):
  U  (S1-S6, auxiliary function, ONE-SIDED):  c_n - main_n <= 3 k q^{n/2}
  L  (S7, twist / Galois cover):              c_n - main_n >= -3 k q^{n/2}
  A  (S8, all n + rationality):               max_i |root_i| <= sqrt(q)
Setup test (S1): can the data be the point count of a genus-1 curve at all?  Needs integer closed-point counts
P_n = (1/n) sum_{d|n} mu(n/d) c_d >= 0 AND Hasse-Deuring |a_1| <= 2 sqrt q (else no curve exists to run S1-S6 on).
A real port must FAIL on Z0 and on FTH, and PASS on CURVE and FTAU."""
import json, math, numpy as np
from fractions import Fraction
def mu(n):
    r, m, d = 1, n, 2
    while d*d <= m:
        if m % d == 0:
            m //= d
            if m % d == 0: return 0
            r = -r
        d += 1
    return -r if m > 1 else r
def ap(p, A, B):
    return p + 1 - (1 + sum(1 + (0 if (x**3 + A*x + B) % p == 0 else (1 if pow((x**3 + A*x + B) % p, (p-1)//2, p) == 1 else -1)) for x in range(p)))
TH, TAU, NMAX = 0.2, 3.0, 40
out = []
for (p, A, B) in [(5, 1, 1), (7, 1, 3), (11, 1, 1), (13, 1, 6)]:
    q = p*p; a1p = ap(p, A, B); aq = a1p*a1p - 2*p          # a_1 over F_{p^2}
    # integer recursions: curve a_n: a_n = aq a_{n-1} - q a_{n-2};  Z0: s_n = q s_{n-1} - q s_{n-2}
    def seq(t, nn):
        s = [2, t]
        for _ in range(nn): s.append(t*s[-1] - q*s[-2])
        return s
    aC, aZ = seq(aq, NMAX), seq(q, NMAX)
    rootsC = np.roots([1, -aq, q]); rootsZ = np.roots([1, -q, q])
    models = {}
    for nm, a, roots in (('CURVE', aC, rootsC), ('Z0', aZ, rootsZ)):
        c = [q**n - a[n] for n in range(NMAX+1)]
        models[nm] = dict(c=[float(c[n]) for n in range(1, NMAX+1)], main=[float(q**n) for n in range(1, NMAX+1)], k=1,
                          roots=[abs(r) for r in roots], c_int=[c[n] for n in range(1, NMAX+1)], a1=a[1], err=[-a[n] for n in range(1, NMAX+1)])
    cC = models['CURVE']['c']
    models['FTAU'] = dict(c=[cC[n-1]*2*math.cos(n*TAU*math.log(q)) for n in range(1, NMAX+1)],
                          main=[q**n*2*math.cos(n*TAU*math.log(q)) for n in range(1, NMAX+1)], k=2,
                          err=[-aC[n]*2*math.cos(n*TAU*math.log(q)) for n in range(1, NMAX+1)], roots=[abs(r) for r in rootsC]*2)
    models['FTH'] = dict(c=[cC[n-1]*(q**(n*TH) + q**(-n*TH)) for n in range(1, NMAX+1)],
                         main=[q**n*(q**(n*TH) + q**(-n*TH)) for n in range(1, NMAX+1)], k=2,
                         err=[-aC[n]*(q**(n*TH) + q**(-n*TH)) for n in range(1, NMAX+1)],
                         roots=[abs(r)*q**TH for r in rootsC] + [abs(r)*q**-TH for r in rootsC])
    for nm, M in models.items():
        e = M['err']; bnd = [3*M['k']*q**((i+1)/2) for i in range(NMAX)]
        Ufail = [i+1 for i in range(NMAX) if e[i] > bnd[i]*(1+1e-12)]; Lfail = [i+1 for i in range(NMAX) if e[i] < -bnd[i]*(1+1e-12)]
        Afail = max(M['roots']) > math.sqrt(q)*(1+1e-12)
        row = dict(p=p, q=q, curve=f'y^2=x^3+{A}x+{B}', model=nm, U_first_fail=Ufail[0] if Ufail else None, U_nfail=len(Ufail),
                   L_first_fail=Lfail[0] if Lfail else None, L_nfail=len(Lfail), A_fail=bool(Afail),
                   max_root_over_sqrtq=max(M['roots'])/math.sqrt(q),
                   max_abs_err_over_bound_n_le_10=max(abs(float(e[i]))/bnd[i] for i in range(10)))
        if 'c_int' in M:
            ci = M['c_int']
            Pn = [Fraction(sum(mu((n)//d)*ci[d-1] for d in range(1, n+1) if n % d == 0), n) for n in range(1, 13)]
            row.update(closed_point_counts_integer_nonneg=all(x.denominator == 1 and x >= 0 for x in Pn),
                       P_1_to_4=[str(x) for x in Pn[:4]], hasse_deuring_curve_possible=abs(M['a1']) <= 2*math.sqrt(q))
        verdict_fail = bool(Ufail or Lfail or Afail)
        row['toy_rejects'] = verdict_fail
        row['expected'] = 'reject' if nm in ('Z0', 'FTH') else 'pass'
        row['correct'] = verdict_fail == (nm in ('Z0', 'FTH'))
        out.append(row)
        print({k: row[k] for k in ('q', 'model', 'U_first_fail', 'L_first_fail', 'A_fail', 'max_root_over_sqrtq', 'toy_rejects', 'correct')},
              {k: row[k] for k in ('closed_point_counts_integer_nonneg', 'P_1_to_4', 'hasse_deuring_curve_possible') if k in row})
summ = dict(all_correct=all(r['correct'] for r in out),
            Z0_rejected_by_U_alone=any(r['U_nfail'] for r in out if r['model'] == 'Z0'),
            Z0_rejected_by_L=all(r['L_first_fail'] for r in out if r['model'] == 'Z0'),
            FTH_rejected_by_U=all(r['U_first_fail'] for r in out if r['model'] == 'FTH'))
print(summ)
json.dump(dict(summary=summ, rows=out), open('falsifier_ge.json', 'w'), indent=1)
