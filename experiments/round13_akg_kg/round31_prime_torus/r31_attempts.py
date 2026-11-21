"""Round 31, task 2 analysis (reads r31_floor.json; recomputes psi for the kernel identity).
 A1 positive-kernel form (PROVED by partial summation; CHECKED here):
      S(L) = sum_{n<=X} Lam(n) n^{-1/2}(L - log n) = int_1^X psi(x) x^{-3/2} (1 + (L - log x)/2) dx,   X = e^L,
    and int_1^X x * x^{-3/2}(1 + (L-log x)/2) dx = 4 sqrt X - L - 4.  Hence
      D(L) = int_1^X (x - psi(x)) x^{-3/2}(1 + (L-log x)/2) dx + [4(e^{-L/2} - 1) + (c1+1) L + C2 - T(L)].
    The floor D >= -B is a ONE-SIDED, positively weighted 'psi(x) <= x + O(sqrt x)' statement.
 A2 Chebyshev-type floors: psi(x) <= c x (c = 1.03883, Rosser-Schoenfeld 1962) gives D >= -4(c-1) e^{L/2} - O(L):
    exponent 1/2 -> through F13W (eps = 1/2) only 'zero-free on Re s > 1' (already known).
 A3 envelope exponents of the window minima (how fast -min D grows) for FTH and EPS, vs. the stretch theta."""
import json, math, numpy as np, mpmath
fl = json.load(open('r31_floor.json')); out = {}
c1, C2 = fl['c1'], fl['C2']
# A1 check at a few L with exact psi (primes up to e^L)
def lam_list(N):
    s = np.ones(N + 1, bool); s[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if s[p]: s[p*p::p] = False
    n, l = [], []
    for p in np.nonzero(s)[0]:
        q = int(p)
        while q <= N: n.append(q); l.append(math.log(p)); q *= int(p)
    o = np.argsort(n); return np.array(n, float)[o], np.array(l)[o]
n, lam = lam_list(int(math.exp(12)) + 1)
rows = []
for L in (4.0, 8.0, 12.0):
    X = math.exp(L); sel = n <= X
    S_direct = float(np.sum(lam[sel]*n[sel]**-0.5*(L - np.log(n[sel]))))
    # integral of a step function psi: on [n_j, n_{j+1}) psi = cumulative; antiderivative of w(x) = x^{-3/2}(1+(L-log x)/2)
    W = lambda x: -x**-0.5*(L - math.log(x))   # d/dx W = w(x)  (checked numerically below)
    pts = np.concatenate([n[sel], [X]]); psi = np.cumsum(lam[sel])
    S_kernel = float(sum(psi[j]*(W(pts[j + 1]) - W(pts[j])) for j in range(len(psi))))
    rows.append(dict(L=L, S_direct=S_direct, S_kernel=S_kernel, diff=abs(S_direct - S_kernel), main=4*math.sqrt(X) - L - 4))
x = mpmath.mpf(7.3); Lm = mpmath.mpf(9)
Wm = lambda x: -x**-0.5*(Lm - mpmath.log(x))
out['A1_antiderivative_check'] = float(mpmath.diff(Wm, x) - x**-1.5*(1 + (Lm - mpmath.log(x))/2))
out['A1_kernel_identity'] = rows; print(rows, out['A1_antiderivative_check'])
# A2
c = 1.03883
out['A2_chebyshev_floor'] = [dict(L=L, floor=-4*(c - 1)*math.exp(L/2) - abs(c1 + 1)*L - 8) for L in (10, 20, 40)]
# A3 envelopes
A3 = {}
for name in ('FTH', 'EPS', 'ZETA', 'DH'):
    Lg = np.array(fl['models'][name]['Lgrid']); Dg = np.array(fl['models'][name]['Dgrid'])
    Ls, env = [], []
    for L0 in np.arange(6, Lg.max()/1.5, 0.25):
        sel = (Lg >= L0) & (Lg <= 1.5*L0)
        m = -Dg[sel].min()
        if m > 0: Ls.append(1.5*L0); env.append(m)
    if len(Ls) > 4:
        sl = np.polyfit(Ls, np.log(env), 1)[0]
        A3[name] = dict(fit_exponent=float(sl), n=len(Ls), last=(Ls[-1], env[-1]))
    else: A3[name] = dict(fit_exponent=None, note='no negative window minimum in range')
out['A3_envelope'] = A3; print(A3)
json.dump(out, open('r31_attempts.json', 'w'), indent=1)
