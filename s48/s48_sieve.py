"""S48: exact-by-sieve psi(x), S(x)=sum_{n<=x} Lambda(n)/sqrt(n), and the
weighted sums needed for G(x), E1(x), log-mean and linear-mean of d_psi.
Segmented numpy sieve; checkpoints at 1e6..1e9 (and dyadic halves x/2).
Convention (ASSUMPTION, Turbo's sign): d_psi(x) = (x - psi(x))/sqrt(x).
 G(x)  = int_1^x (psi(t)-t) t^{-3/2} dt = 2 S(x) - 2 psi(x)/sqrt(x) - 2 sqrt(x) + 2
 E1(x) = -d_psi(x) - (1/2) int_1^x d_psi dt/t = S(x) - 2 sqrt(x) + 1   (exact)
 Lin(x)= int_1^x d_psi(t) dt = (2/3)(x^{3/2}-1) - 2 sqrt(x) psi(x) + 2 R(x),  R=sum Lambda(n) sqrt(n)
Label: NUMERICAL (float64 sums; checked against mpmath at 1e6)."""
import numpy as np, math, json, sys
XMAX = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**9
cks = sorted(set([10**k for k in range(6, 11) if 10**k <= XMAX] + [10**k // 2 for k in range(6, 11) if 10**k <= XMAX]))
r = int(math.isqrt(XMAX)) + 1
base = np.ones(r + 1, bool); base[:2] = False
for i in range(2, int(r**0.5) + 1):
    if base[i]: base[i*i::i] = False
bp = np.nonzero(base)[0]
# prime powers p^k, k>=2 (few): collect
pp = []
for p in bp:
    q = int(p) * int(p)
    while q <= XMAX:
        pp.append((q, math.log(int(p)))); q *= int(p)
pp.sort()
ppn = np.array([a for a, _ in pp], dtype=np.int64); ppl = np.array([b for _, b in pp])
psi = 0.0; S = 0.0; R = 0.0
cS = 0.0  # Kahan compensation for S
out = {}
SEG = 1 << 24
lo = 2; ci = 0
def add(vals_n, vals_l):
    global psi, S, R
    psi += math.fsum(vals_l); S += math.fsum(vals_l / np.sqrt(vals_n)); R += math.fsum(vals_l * np.sqrt(vals_n))
while lo <= XMAX and ci < len(cks):
    hi = min(lo + SEG, XMAX + 1)
    seg = np.ones(hi - lo, bool)
    for p in bp:
        p = int(p)
        if p * p >= hi: break
        st = max(p * p, ((lo + p - 1) // p) * p)
        seg[st - lo::p] = False
    primes = np.nonzero(seg)[0].astype(np.int64) + lo
    m = (ppn >= lo) & (ppn < hi)
    alln = np.concatenate([primes, ppn[m]]); alll = np.concatenate([np.log(primes.astype(float)), ppl[m]])
    order = np.argsort(alln); alln = alln[order]; alll = alll[order]
    # handle checkpoints inside segment
    while ci < len(cks) and cks[ci] < hi:
        x = cks[ci]; k = np.searchsorted(alln, x, side='right')
        add(alln[:k].astype(float), alll[:k]); alln = alln[k:]; alll = alll[k:]
        sx = math.sqrt(x)
        G = 2*S - 2*psi/sx - 2*sx + 2
        dpsi = (x - psi)/sx
        E1 = S - 2*sx + 1
        Lin = (2/3)*(x**1.5 - 1) - 2*sx*psi + 2*R
        out[x] = dict(x=x, psi=psi, S=S, d_psi=dpsi, G=G, E1=E1, logmean=-G/math.log(x), linmean=Lin/x)
        print(json.dumps(out[x]), flush=True); ci += 1
    add(alln.astype(float), alll)
    lo = hi
json.dump(out, open('sieve_%d.json' % XMAX, 'w'), indent=1)
