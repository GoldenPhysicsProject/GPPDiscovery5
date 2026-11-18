"""Round 13h test: Connes' cutoff = Weil form restricted to test functions supported in [lambda^-1/2, lambda^1/2]
(our width L = log lambda). lam_min(W_L) for zeta, Epstein (own arch) and spliced Beurling B'' (zeta's arch borrowed;
real primes < 55, Beurling primes >= 55, quartet 0.6+-28i) at L = 4 ... 8. For L <= log 55 B'' == zeta exactly."""
import numpy as np, json, math
import harness as h
Z = h.zeta_model(); E = h.epstein_model(); B = h.beurling_model(x0=55.0)
m = Z['n'] < 55
S = dict(Z, name='SPLICE', n=np.concatenate([Z['n'][m], B['n']]), c=np.concatenate([Z['c'][m], B['c']]))
o = np.argsort(S['n']); S['n'] = S['n'][o]; S['c'] = S['c'][o]
out = {}
for L in (4.0, 4.5, 5.0, 6.0, 7.0, 8.0):
    for M in (int(10*L), int(14*L)):
        row = {}
        for mo in (Z, E, S):
            a, p, pr, _ = h.operators(mo, L, M); ev = np.linalg.eigvalsh(a+p-pr)
            row[mo['name']] = dict(lam_min=float(ev[0]), n_neg=int(np.sum(ev < -1e-8)),
                                   signed_bound=float(np.linalg.eigvalsh(pr-p)[-1]))
        out[f'L{L}_M{M}'] = row; print(L, M, json.dumps(row), flush=True)
json.dump(out, open('connes_cutoff.json', 'w'), indent=1)
