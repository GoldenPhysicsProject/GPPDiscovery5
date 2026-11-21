"""Round 34 analysis of the ratio-class layers (r34_layers.c output) for the real mu, the Codex sparse twist, the
Selberg-Delange main-term shape of E_X(z), the prime-class Dirichlet-series main term, and the twist's identity ledger."""
import json, math, numpy as np, mpmath
Xs = [10000, 30000, 100000, 300000, 1000000, 3000000, 10000000]
c = 1/-1.4603545088095868
def load(X, t): d = json.load(open(f'{"lt" if t else "l"}_{X}.json')); return d, np.array([e[0] for e in d['Ek']]), max(abs(e[1]) for e in d['Ek'])
# Selberg-Delange constants A(z) = prod_p (1 - z/p)(1 - 1/p)^{-z},  J(z) = prod_p (1 - 2z/p + 1/p^2)/(1 - z/p)^2
N = 10**7; sv = np.ones(N + 1, bool); sv[:2] = False
for p in range(2, int(N**0.5) + 1):
    if sv[p]: sv[p*p::p] = False
pr = np.nonzero(sv)[0].astype(float)
def A(z): return float(np.exp(np.sum(np.log1p(-z/pr) - z*np.log1p(-1/pr))))
def J(z): return float(np.exp(np.sum(np.log1p(-2*z/pr + 1/pr**2) - 2*np.log1p(-z/pr))))
def rgam(z): return float(mpmath.rgamma(-z))
res = dict(rows=[], twist_rows=[], curve={}, prime_class={}, ledger={})
for X in Xs:
    for t in (0, 1):
        d, E, imag = load(X, t); L = d['L']; Eabs = np.abs(E).sum()
        row = dict(X=X, L=L, E0_H=E[0], E1_prime=E[1], E2=E[2], E3=E[3], sum_k_ge2=float(E[2:].sum()), R=float(E[1:].sum()),
                   G_E1=float(E.sum()), G_direct=d['E_at_1'], E_minus1=d['E_at_minus1'], kmax=int(np.argmax(np.abs(E))),
                   maxabs=float(np.abs(E).max()), sum_abs=float(Eabs), cancellation=float(E.sum()/Eabs), dft_imag_max=imag,
                   E1_over_sqrtX_logX=float(abs(E[1])/(math.sqrt(X)/L)), maxabs_over_X=float(np.abs(E).max()/X),
                   Ek=E[:14].tolist())
        (res['twist_rows'] if t else res['rows']).append(row)
    # E_X(z) on [-1,1] and the Selberg-Delange main-term shape  L * 4X A^2 J (log X)^{-2z-2} / Gamma(-z)^2
    d, E, _ = load(X, 0); L = d['L']; zs = np.linspace(-1, 1, 81)
    Ez = np.array([np.polyval(E[::-1], z) for z in zs])
    main = np.array([L*4*X*A(z)**2*J(z)*L**(-2*z - 2)*rgam(z)**2 for z in zs])
    res['curve'][X] = dict(z=zs.tolist(), E=Ez.tolist(), SD_main=main.tolist(), min_E_on_interval=float(Ez.min()))
# prime class: exact (from layers) vs the Dirichlet-series main term of S_p(y) = sum_{g<=y, p!|g} mu(g)^2/g
#   S_q(y) = (6/pi^2) prod_{p|q} p/(p+1) * [log y + gamma - 2 zeta'(2)/zeta(2) + sum_{p|q} log p/(p+1)] + O(y^{-1/2})
C0 = float(mpmath.euler - 2*mpmath.zeta(2, derivative=1)/mpmath.zeta(2))
for X in (10000, 100000, 1000000, 10000000):
    L = math.log(X); p = pr[pr <= X]; y = X//p
    Smain = 6/math.pi**2*p/(p + 1)*(np.log(y) + C0 + np.log(p)/(p + 1))
    E1main = float(-2*np.sum((L - np.log(p))/np.sqrt(p)*(Smain - c)))
    E1 = [r for r in res['rows'] if r['X'] == X][0]['E1_prime']
    asym = -2*math.sqrt(X)/L*(6/math.pi**2*(16 + 4*C0) - 4*c)
    res['prime_class'][X] = dict(exact=E1, S_main_term_model=E1main, leading_asymptotic=asym, ratio_exact_to_leading=E1/asym)
# S_q main term spot check
mu = np.ones(200001); sq = np.ones(200001, bool)
for q in pr[pr <= 450]: sq[int(q*q)::int(q*q)] = False
chk = {}
for q in (2, 3, 30, 210):
    for Y in (1000, 100000):
        g = np.arange(1, Y + 1); m = sq[1:Y + 1] & (np.gcd(g, q) == 1); ex = float(np.sum(1/g[m]))
        ps = [pp for pp in (2, 3, 5, 7) if q % pp == 0]; fac = np.prod([pp/(pp + 1) for pp in ps])
        mt = 6/math.pi**2*fac*(math.log(Y) + C0 + sum(math.log(pp)/(pp + 1) for pp in ps))
        chk[f'q{q}_Y{Y}'] = dict(exact=ex, main=mt, diff=ex - mt, diff_times_sqrtY=(ex - mt)*math.sqrt(Y))
res['S_main_check'] = chk; res['C0'] = C0
# twist ledger at X = 1e6: which identities the sparse flip violates
X = 10**6; mu = np.ones(X + 1, np.int8); mu[0] = 0; s2 = np.ones(X + 1, bool); s2[:2] = False
for p in range(2, X + 1):
    if s2[p]:
        s2[p*p::p] = False; mu[p::p] *= -1; mu[p*p::p*p] = 0
isp = s2; flip = np.zeros(X + 1, bool); flip[X//2 + 1: 3*X//4 + 1] = isp[X//2 + 1: 3*X//4 + 1]
mut = mu.copy(); mut[flip] *= -1; n = np.arange(1, X + 1)
for nm, m_ in (('mu', mu), ('twist', mut)):
    a = m_[1:]/np.sqrt(n)
    res['ledger'][nm] = dict(sum_mu_over_n=float(np.sum(m_[1:]/n)), meissel_sum_mu_floor=int(np.sum(m_[1:].astype(np.int64)*(X//n))),
                             F1_centered=float(a.sum() - c), M_X=int(m_[1:].sum()), n_flipped=int(flip.sum()),
                             max_divisor_identity_violation=2 if nm == 'twist' else 0)
res['ledger']['shift_2sum_p^-1/2'] = float(2*np.sum(1/np.sqrt(n[flip[1:]])))
res['ledger']['shift_2sum_1/p'] = float(2*np.sum(1/n[flip[1:]]))
json.dump(res, open('r34_results.json', 'w'), indent=1, default=float)
for r in res['rows'] + res['twist_rows']:
    print('X=%8d E0=%8.2f E1=%10.2f E2=%12.4g R=%10.2f G=%9.4f E(-1)=%.4g kmax=%d max=%.4g max/X=%.3f canc=%.2e |E1|/(sqrtX/L)=%.1f' %
          (r['X'], r['E0_H'], r['E1_prime'], r['E2'], r['R'], r['G_E1'], r['E_minus1'], r['kmax'], r['maxabs'], r['maxabs_over_X'], r['cancellation'], r['E1_over_sqrtX_logX']))
print(json.dumps(res['prime_class'], indent=0)); print(json.dumps(res['ledger'], indent=0)); print(json.dumps(res['S_main_check'], indent=0))
for X in Xs: cv = res['curve'][X]; print(X, 'min E on [-1,1] = %.4g' % cv['min_E_on_interval'])
