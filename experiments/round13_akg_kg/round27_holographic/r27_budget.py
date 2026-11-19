"""Round 27 conservation-law hunt, hardest check: does the per-prime charge (Lambda(p^k) = log p exactly, i.e.
'each prime counted once') carry the positivity, or is positivity a global balance?
(i)  'no budget': scale the prime side, Lambda -> (1+eps) Lambda, and find the eps-window in which the sector Weil
     form stays PSD (zeta, L = 1, 1.5, where lam_min is above the double floor).
(ii) truncation: keep only primes p <= P (all their powers) inside the window; is the truncated form PSD?
(iii) floor: prime side vs zero side (300 zeros) for zeta at L = 2, 5 after the orbit fix (sets the numerical floor).
NUMERICAL. No zeros on the prime side."""
import numpy as np, json, math
from r22_models import lam_table
from r27_sectors import gram_prime, gram_zero, klein
out = {}
lam = lam_table('ZETA', 400)
for L in (1.0, 1.5):
    M = int(3 * L + 6)
    for par, nm in ((1, 'even'), (-1, 'odd')):
        W1 = gram_prime('ZETA', L, M, par, lam=lam)
        W0 = gram_prime('ZETA', L, M, par, lam=np.zeros_like(lam))      # pole + arch only
        P = W0 - W1                                                    # prime part (subtracted)
        sc = np.abs(W1).max()
        lm = lambda e: float(np.linalg.eigvalsh(W0 - (1 + e) * P)[0] / sc)
        eps = np.concatenate([-np.logspace(-6, 0, 61)[::-1], [0.0], np.logspace(-6, 0, 61)])
        vals = [lm(e) for e in eps]
        ok = [e for e, v in zip(eps, vals) if v > 0]
        # truncation: primes p <= Pmax only
        trunc = {}
        for Pm in (0, 2, 3, 5, 7, 11, 13, 17, 19):
            lt = lam.copy()
            for n in range(2, len(lt)):
                if lt[n] > 0:
                    p = round(math.exp(lt[n]))
                    if p > Pm: lt[n] = 0.0
            Wt = gram_prime('ZETA', L, M, par, lam=lt)
            trunc[str(Pm)] = float(np.linalg.eigvalsh(Wt)[0] / sc)
        out[f'L{L}_{nm}'] = dict(lam_min_rel=lm(0.0), eps_window=[float(min(ok)), float(max(ok))] if ok else None,
                                 lam_min_vs_eps={f'{e:.3g}': v for e, v in zip(eps, vals) if abs(e) in (1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1) or e == 0},
                                 truncation_lam_min_rel=trunc, primes_in_window=[p for p in (2, 3, 5, 7, 11, 13, 17, 19) if math.log(p) < L])
        print(L, nm, json.dumps(out[f'L{L}_{nm}']), flush=True)
zz = json.load(open('zeta_zeros_300.json'))
for L in (2.0, 5.0):
    M = int(3 * L + 6)
    for par, nm in ((1, 'even'), (-1, 'odd')):
        Wp = gram_prime('ZETA', L, M, par); Wz = gram_zero(klein([complex(t) for t in zz]), L, M, par)
        out[f'floor_L{L}_{nm}'] = dict(maxdiff_rel=float(np.abs(Wp - Wz).max() / np.abs(Wz).max()),
                                       opnorm_diff_rel=float(np.linalg.norm(Wp - Wz, 2) / np.abs(Wz).max()),
                                       lam_min_prime_rel=float(np.linalg.eigvalsh(Wp)[0] / np.abs(Wz).max()),
                                       lam_min_zero_rel=float(np.linalg.eigvalsh(Wz)[0] / np.abs(Wz).max()))
        print('floor', L, nm, out[f'floor_L{L}_{nm}'], flush=True)
json.dump(out, open('r27_budget.json', 'w'), indent=1)
