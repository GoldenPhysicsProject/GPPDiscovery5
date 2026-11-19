"""Task 1: S-fraction / Cauer ladder for each model at two precisions; report stable digits and first non-positive index."""
import json, sys, time
import mpmath as mp
from r20_core import *
K = int(sys.argv[1]) if len(sys.argv) > 1 else 70
D1, D2 = (int(x) for x in (sys.argv[2:4] if len(sys.argv) > 3 else (400, 700)))
models = sys.argv[4].split(',') if len(sys.argv) > 4 else ['ZETA', 'FTH', 'DH', 'EPS']
out = {}
for name in models:
    res = {}
    for D in (D1, D2):
        t = time.time()
        f = F_coeffs(name, K, D)
        mp.mp.dps = D
        mu = m_series(f, K)
        c0, a = s_fraction(mu, K - 1)
        G, C = ladder(c0, a)
        res[D] = dict(f=f, mu=mu, c0=c0, a=a, G=G, C=C)
        print(name, D, 'done', round(time.time() - t, 1), 's', flush=True)
    A1, A2 = res[D1]['a'], res[D2]['a']
    mp.mp.dps = D2
    digits = [float(-mp.log10(abs(x - y)/abs(y))) if x != y else float(D1) for x, y in zip(A1, A2)]
    stable = [k for k, d in enumerate(digits) if d >= 10]
    nstable = 0
    while nstable < len(digits) and digits[nstable] >= 10: nstable += 1
    a = A2[:nstable]
    neg = next((k + 1 for k, x in enumerate(a) if x <= 0), None)
    row = dict(c0=mp.nstr(res[D2]['c0'], 30), n_stable_a=nstable, first_nonpositive_a_index=neg,
               a=[mp.nstr(x, 20) for x in a], stable_digits=[round(d, 1) for d in digits],
               G=[mp.nstr(x, 15) for x in res[D2]['G'][:nstable//2 + 1]],
               C=[mp.nstr(x, 15) for x in res[D2]['C'][:nstable//2]],
               moments_c=[mp.nstr(((-1)**j)*res[D2]['mu'][j], 20) for j in range(12)])
    out[name] = row
    print(name, 'stable', nstable, 'first<=0 at', neg, 'a[:12]=', [mp.nstr(x, 6) for x in a[:12]], flush=True)
    json.dump(out, open(f'r20_ladder_{"_".join(models)}.json', 'w'), indent=1)
