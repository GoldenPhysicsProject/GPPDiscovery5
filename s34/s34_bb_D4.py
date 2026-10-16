"""Grok S31 H31a: rigorous upper bound on A = sup_torus (B - m) for the S30 D=4 Fejér–Riesz torus lift.
Reuses s28_torus_bb.bb (float64 B&B, 2nd-order Taylor/Lagrange remainder, explicit a-priori rounding margin).
Coefficients a are rounded to 10-significant-digit decimal strings a_q (exactly what an Arb certificate reads).
Invisibility of every added character checked in Arb (200 bits): |k.log p| > 2L strictly.
Prime set = all primes p with log p < 2L, so the torus covers the whole visible comb (no non-torus remainder).
Usage: python s31_bb_D4.py KEY [margin]   (KEY e.g. L1.0_D4).  RH OPEN."""
import sys, os, json, math, time
import numpy as np
from flint import arb, ctx
ctx.prec = 200
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's28'))
from s28_torus_bb import bb
key = sys.argv[1]; margin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
SRC = os.environ.get('SRC', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's30', 'h30a_fejer_riesz.json'))
src = json.load(open(SRC))[key]
L = src['L']; ps = src['primes']; n = len(ps)
# prime set completeness
allp = [q for q in range(2, 200) if all(q % r for r in range(2, int(q**.5) + 1)) and float(arb(q).log()) < 2 * L]
if os.environ.get('SUBSET') != '1': assert allp == ps, (allp, ps)
else: print('SUBSET mode: remaining primes', [q for q in allp if q not in ps], 'must be bounded separately (s31_merge_p11.py)', flush=True)
# exact visible coefficients and invisibility check
lps = [arb(p).log() for p in ps]; twoL = 2 * arb(str(L))
for (i, k, c) in src['vis']:
    assert k * lps[i] < twoL  # arb: certainly less
    assert abs(arb(c) - 2 * lps[i] / arb(ps[i]) ** (arb(k) / 2)) < arb('1e-14')
minsep = None
for kv in src['K']:
    u = abs(sum((int(kk) * lp for kk, lp in zip(kv, lps)), arb(0)))
    assert u > twoL, ('VISIBLE', kv); minsep = float((u - twoL).lower()) if minsep is None or float((u - twoL).lower()) < minsep else minsep
print(f'{key}: L={L} primes={ps} chars={len(src["K"])} all invisible; min(|k.logp|-2L)={minsep:.3e}', flush=True)
rows, coefs = [], []
for (i, k, c) in src['vis']:
    kv = [0] * n; kv[i] = k; rows.append(kv); coefs.append(float((2 * lps[i] / arb(ps[i]) ** (arb(k) / 2)).mid()))
a_q = ['%.9e' % x for x in src['a_lp']]
for kv, x in zip(src['K'], a_q): rows.append(kv); coefs.append(-float(x))
Kmat = np.array(rows, float); coef = np.array(coefs)
A_est = src['A_true']; target = float('%.4f' % (A_est + margin)); t0 = time.time()
maxb = int(os.environ.get('MAXBOXES', '12000000'))
while True:
    r, best, nproc, Er = bb(Kmat, coef, n, target, maxboxes=maxb)
    if r is not None: break
    if best > target: target = float('%.4f' % (best + margin))
    else: target = float('%.4f' % (target + margin))
    print('raise target', target, flush=True)
print(f'{key}: sup_torus (B-m) <= {target} (best sample {best:.6f}); boxes {nproc}; Eround {Er:.2e} [{time.time()-t0:.0f}s]', flush=True)
out = dict(L=L, primes=ps, vis=[list(v) for v in src['vis']], K=src['K'], a_q=a_q, A_cert=repr(target), A_best_sample=best,
           per_prime=src['per_prime'], lam_max_floor=src['lam_max'], min_invisible_gap=str(minsep),
           method='float64 B&B (s28_torus_bb), 2nd-order Taylor, explicit rounding margin; libm cos/sin within 2 ulp assumed',
           boxes=nproc, Eround=Er, elapsed_s=time.time() - t0)
json.dump(out, open(f'torus_{key}_cert.json', 'w'), indent=1)
