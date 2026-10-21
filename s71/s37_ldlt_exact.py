"""Grok S37 H37b: SECOND, INDEPENDENT positive-definiteness certificate for a dumped Gram midpoint matrix Atil
(H3 of the S35 conditional theorem).  Shares NO code with s36_cholres.py and uses NO Arb / interval arithmetic:
everything is exact integer (Python int + python-flint fmpz_mat) and exact rationals (fractions.Fraction).

Method (exact-integer fixed-point LDL^T + exact residual + Frobenius/Weyl):
  0. sha256 of the dump is checked against the expected value (argument).
  1. Each decimal string x_ij is parsed EXACTLY (integer mantissa, integer exponent) and rounded to
     N_ij = round(x_ij * 2^S).  |x_ij - N_ij/2^S| <= 2^-(S+1).
     Dump-string error: the strings are arb .str(radius=False) of the 384-bit midpoints, printed with
     nd significant digits; we take |Atil_ij - x_ij| <= 1 unit in the last printed digit = 10^(e_last)
     (conservative: correct rounding would give 1/2 unit).
  2. c0 = floor(l0 * 2^S).  Fixed-point (non-rigorous, any result is fine) LDL^T of W = sym(N) - c0 I
     with L unit lower triangular (scale 2^S) and D diagonal (scale 2^S), integer floor arithmetic.
     If some pivot d_k <= 0 the attempt FAILS (no claim).
  3. RIGOROUS, EXACT: E = 2^(2S) (N - c0 I) - L_int diag(D_int) L_int^T  (fmpz_mat, exact integers).
     All D_int > 0 (checked exactly) => L D L^T >= 0, hence for all real x:
       x^T (N/2^S - c0/2^S I) x >= -||E||_2 / 2^(3S) |x|^2 >= -||E||_F / 2^(3S) |x|^2   (E need not be symmetric).
  4. Weyl:  lambda_min(sym Atil) >= c0/2^S - ||E||_F/2^(3S) - ||Atil - N/2^S||_F,
     ||Atil - N/2^S||_F <= sqrt(sum_ij (10^(e_last,ij) + 2^-(S+1))^2)  (computed as an exact rational, sqrt rounded up).
  The returned bound is an exact rational rounded DOWN to 6 significant digits.
Status: computer-assisted, second implementation, NOT Lean.  RH OPEN.
usage: s37_ldlt_exact.py dump.json sha256 S l0[,l0b,...]"""
import sys, json, time, hashlib, math
from fractions import Fraction
import numpy as np
from flint import fmpz_mat

fn, sha_exp, S = sys.argv[1], sys.argv[2], int(sys.argv[3]); l0s = sys.argv[4].split(',')
t0 = time.time()
raw = open(fn, 'rb').read(); sha = hashlib.sha256(raw).hexdigest()
assert sha == sha_exp, ('SHA256 MISMATCH', sha)
rows = json.loads(raw); K = len(rows); assert all(len(r) == K for r in rows)

def parse(s):
    """exact: returns (mantissa int, exponent int) with value = m*10^e, and e_last = exponent of last printed digit."""
    s = s.strip(); neg = s.startswith('-'); s = s.lstrip('+-')
    mant, _, ex = s.lower().partition('e'); ex = int(ex) if ex else 0
    ip, _, fp = mant.partition('.')
    m = int((ip + fp) or '0'); e = ex - len(fp)
    return (-m if neg else m), e

def round_scaled(m, e, S):  # round(m*10^e*2^S) exactly, ties away (any rounding with |err|<=1/2 is fine)
    if e >= 0: return m * 10 ** e << S
    num = m << S; den = 10 ** (-e)
    q, r = divmod(abs(num), den); q += (2 * r >= den)
    return q if num >= 0 else -q

N = [[0] * K for _ in range(K)]; elast = []
for i in range(K):
    for j in range(K):
        m, e = parse(rows[i][j]); N[i][j] = round_scaled(m, e, S); elast.append(e)
del rows, raw
# Frobenius bound on ||Atil - N/2^S||: entries <= 10^e + 2^-(S+1); group by e
from collections import Counter
cnt = Counter(elast); half = Fraction(1, 2 ** (S + 1))
fro2 = sum(c * (Fraction(10) ** e + half) ** 2 for e, c in cnt.items())
def sqrt_up(q):  # rational upper bound of sqrt(q), q Fraction >= 0
    sc = 2 ** 2000; v = math.isqrt(q.numerator * sc * sc // q.denominator + 1) + 1
    return Fraction(v, sc)
dpert = sqrt_up(fro2)
print(f'K={K} S={S} sha256 OK; max e_last={max(cnt)}; ||Atil-N/2^S||_F <= {float(dpert):.3e} [{time.time()-t0:.0f}s]', flush=True)

Nobj = np.empty((K, K), dtype=object)
for i in range(K):
    for j in range(K): Nobj[i, j] = N[i][j]
res = {'dump': fn, 'sha256': sha, 'K': K, 'S': S, 'pert_F': float(dpert), 'tries': []}
for l0s_ in l0s:
    l0 = Fraction(l0s_); c0 = math.floor(l0 * 2 ** S); t1 = time.time()
    W = (Nobj + Nobj.T) // 2
    for k in range(K): W[k, k] -= c0
    Lm = np.zeros((K, K), dtype=object); Lm[:] = 0; D = [0] * K; fail = None; one = 1 << S
    for k in range(K):
        d = W[k, k]
        if d <= 0: fail = k; break
        D[k] = d; Lm[k, k] = one
        if k + 1 < K:
            col = W[k + 1:, k]
            l = (col << S) // d
            Lm[k + 1:, k] = l
            W[k + 1:, k + 1:] -= np.multiply.outer(l, W[k, k + 1:]) >> S
    if fail is not None:
        print(f' l0={l0s_}: fixed-point LDL^T FAILED at pivot {fail} (W - l0 I numerically indefinite) [{time.time()-t1:.0f}s]', flush=True)
        res['tries'].append({'l0': l0s_, 'ldlt': 'fail', 'pivot': fail}); continue
    assert all(int(x) > 0 for x in D)
    minD = min(D)
    X = Lm * np.array(D, dtype=object)[None, :]
    t2 = time.time()
    P = fmpz_mat([[int(x) for x in r] for r in X.tolist()]) * fmpz_mat([[int(x) for x in r] for r in Lm.T.tolist()])
    sh = 2 * S; ssq = 0
    for i in range(K):
        for j in range(K):
            v = (N[i][j] - (c0 if i == j else 0)) << sh
            v -= int(P[i, j]); ssq += v * v
    normE = Fraction(math.isqrt(ssq) + 1, 2 ** (3 * S))
    lam = Fraction(c0, 2 ** S) - normE - dpert
    ok = lam > 0
    # round down to 6 significant digits
    ex = math.floor(math.log10(float(lam))) if ok else 0
    lam_dn = Fraction(math.floor(lam / Fraction(10) ** (ex - 5)), 1) * Fraction(10) ** (ex - 5) if ok else None
    print(f' l0={l0s_}: minD={float(Fraction(minD, 2**S)):.3e} ||E||_F<={float(normE):.3e} => lambda_min(sym Atil) >= '
          f'{float(lam):.6e} {"OK" if ok else "NO"} [LDLt {t2-t1:.0f}s, residual {time.time()-t2:.0f}s]', flush=True)
    res['tries'].append({'l0': l0s_, 'minD': float(Fraction(minD, 2 ** S)), 'normE_F': float(normE), 'lam_lower': float(lam),
                         'lam_lower_rounded_down': (str(lam_dn.numerator) + '/' + str(lam_dn.denominator)) if ok else None, 'ok': bool(ok)})
    if ok and 'lam_lower' not in res: res['lam_lower'] = float(lam_dn); res['l0'] = l0s_
out = fn.split('/')[-1].replace('.json', f'_ldltexact_S{S}.json')
json.dump(res, open(sys.argv[5] if len(sys.argv) > 5 else out, 'w'), indent=1)
print('DONE', time.time() - t0, flush=True)
