"""Round 29 (C): multiplicity-one integer test. Integer N = C(2n,n) (Chebyshev). Forced divisibility: every prime in (n,2n]
divides N once, so ASK := theta(2n) - theta(n). Size: MAIN := log C(2n,n) <= 2n log 2. The full Legendre sum equals MAIN
exactly (unique factorisation). The slack MAIN - ASK grows like ~0.39 n, never like sqrt(n): at multiplicity one the
trivial loss is of size x. CHECKED (exact integer/float arithmetic). Run: python3 chebyshev_slack.py"""
import math, numpy as np
def primes(n):
    s = np.ones(n + 1, bool); s[:2] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = 0
    return np.nonzero(s)[0]
for n in (10**3, 10**4, 10**5, 10**6):
    P = primes(2 * n); th = lambda x: float(np.sum(np.log(P[P <= x])))
    ask = th(2 * n) - th(n); main = math.lgamma(2 * n + 1) - 2 * math.lgamma(n + 1)
    def vp(p):
        s, k = 0, p
        while k <= 2 * n: s += (2 * n) // k - 2 * (n // k); k *= p
        return s
    full = sum(math.log(p) * vp(int(p)) for p in P)
    print(f"n={n:>8}  ASK={ask:12.1f}  MAIN={main:12.1f}  full Legendre={full:12.1f}  "
          f"slack={main-ask:10.1f}  slack/n={(main-ask)/n:.3f}  slack/sqrt(2n)={(main-ask)/math.sqrt(2*n):7.2f}")
