"""Grok S33: merge a 4D torus cert (primes 2,3,5,7; s31_bb_D4) with S28 per-prime lifts for extra primes (11, 13, ...).
B - m is separable over prime groups: sup(B-m) = sup_{T^4}(B4-m4) + sum_p sup_theta (sum_k c_pk cos k theta - m_p(theta)).
Each kappa_p bounded in Arb on a 40000-interval grid of [0,pi] (as s31_merge_p11 / s28_onestroke_fast --lift).
Invisibility of every added character Arb-checked at the target L; the visible ks are checked visible.
Usage: python s33_merge_pp.py torus4.json L out.json p1:liftjson1 [p2:liftjson2 ...]   RH OPEN."""
import json, sys, math
from flint import arb, ctx
ctx.prec = 128
t4 = json.load(open(sys.argv[1])); L = arb(sys.argv[2]); out = sys.argv[3]
lps4 = [arb(p).log() for p in t4['primes']]
for kv in t4['K']:
    assert abs(sum((int(k) * lp for k, lp in zip(kv, lps4)), arb(0))) > 2 * L, ('visible at target L', kv)
extra = [(int(s.split(':')[0]), s.split(':')[1]) for s in sys.argv[4:]]
n0 = len(t4['primes']); ntot = n0 + len(extra)
K = [list(k) + [0] * len(extra) for k in t4['K']]; a_q = list(t4['a_q']); kaps = {}
A = arb(t4['A_cert'])
for j, (p, lf) in enumerate(extra):
    d = json.load(open(lf))['primes'][str(p)]; lp = arb(p).log()
    ks = d['ks']; cs = [2 * lp / arb(p) ** (arb(k) / 2) for k in ks]
    inv = d['inv']; aq = ['%.12e' % x for x in d['a']]; av = [arb(x) for x in aq]
    for k in inv: assert arb(k) * lp > 2 * L
    for k in ks: assert arb(k) * lp < 2 * L
    ng = 40000; hw = arb.pi() / (2 * ng); kap = None
    for i in range(ng):
        th = arb(((2 * i + 1) * hw).mid(), (hw * 1.0001).mid())
        v = sum((c * (k * th).cos() for k, c in zip(ks, cs)), arb(0)) - sum((x * (k * th).cos() for k, x in zip(inv, av)), arb(0))
        u = v.upper(); kap = u if kap is None or u > kap else kap
    kaps[p] = float(kap); A += arb(kap.upper())
    for k, x in zip(inv, aq):
        kv = [0] * ntot; kv[n0 + j] = k; K.append(kv); a_q.append(x)
Af = math.ceil(float(A.upper()) * 1e4 + 1e-9) / 1e4
res = dict(L=sys.argv[2], primes=t4['primes'] + [p for p, _ in extra], K=K, a_q=a_q, A_cert=repr(Af), A4_cert=t4['A_cert'],
           kappa_upper={str(k): v for k, v in kaps.items()}, source4=sys.argv[1], sources={str(p): f for p, f in extra},
           method='4D torus B&B (s31_bb_D4) + Arb-grid per-prime kappa_p (separable sum)')
json.dump(res, open(out, 'w'), indent=1)
print(f'L={sys.argv[2]} A4={t4["A_cert"]} kappas={kaps} A_cert={Af}')
