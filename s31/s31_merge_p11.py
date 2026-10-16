"""Grok S31: merge a 4D torus cert (primes 2,3,5,7) with the S28 per-prime lift for p=11 at L=1.25.
B - m is separable: sup_{T^5}(B-m) = sup_{T^4}(B4-m4) + sup_theta (c11 cos theta - m11(theta)).
kappa_11 bounded in Arb on a 40000-interval grid of [0,pi] (same as s28_onestroke_fast --lift).
Writes a --torus JSON for s28_onestroke_fast.py with primes [2,3,5,7,11] and A_cert = A4 + kappa11 (upper, 4 decimals up). RH OPEN."""
import json, sys, math
from flint import arb, ctx
ctx.prec = 128
t4 = json.load(open(sys.argv[1])); lj = json.load(open(sys.argv[2])); L = arb(sys.argv[3]); out = sys.argv[4]
d = lj['primes']['11']; p = 11; lp = arb(p).log()
ks = d['ks']; cs = [2 * lp / arb(p) ** (arb(k) / 2) for k in ks]
inv = d['inv']; a_q = ['%.12e' % x for x in d['a']]; av = [arb(x) for x in a_q]
for k in inv: assert arb(k) * lp > 2 * L
for k in ks: assert arb(k) * lp < 2 * L
ng = 40000; hw = arb.pi() / (2 * ng); kap = None
for i in range(ng):
    th = arb(((2 * i + 1) * hw).mid(), (hw * 1.0001).mid())
    v = sum((c * (k * th).cos() for k, c in zip(ks, cs)), arb(0)) - sum((x * (k * th).cos() for k, x in zip(inv, av)), arb(0))
    u = v.upper(); kap = u if kap is None or u > kap else kap
kap = float(kap)
A4 = float(t4['A_cert']); A = math.ceil((A4 + kap) * 1e4 + 1e-9) / 1e4
K = [list(k) + [0] for k in t4['K']] + [[0, 0, 0, 0, k] for k in inv]
res = dict(L=str(L.mid()), primes=t4['primes'] + [11], K=K, a_q=list(t4['a_q']) + a_q, A_cert=repr(A), A4_cert=t4['A_cert'],
           kappa11_upper=kap, source4=sys.argv[1], method='4D torus B&B (s31_bb_D4) + Arb-grid kappa_11 (separable sum)')
json.dump(res, open(out, 'w'), indent=1)
print(f'A4={A4} kappa11<={kap:.6f} A_cert={A}')
