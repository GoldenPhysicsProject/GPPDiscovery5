"""Grok S33: merge two disjoint-prime torus certs (separable sum): A_cert = A_1 + A_2 (each from s31_bb_D4 B&B).
Usage: python s33_merge_torus2.py t1.json t2.json L out.json.  Invisibility of every character re-checked in Arb at L.  RH OPEN."""
import json, sys, math
from flint import arb, ctx
ctx.prec = 128
t1 = json.load(open(sys.argv[1])); t2 = json.load(open(sys.argv[2])); L = arb(sys.argv[3])
ps = t1['primes'] + t2['primes']; lps = [arb(p).log() for p in ps]; n1, n2 = len(t1['primes']), len(t2['primes'])
K = [list(k) + [0] * n2 for k in t1['K']] + [[0] * n1 + list(k) for k in t2['K']]
for kv in K: assert abs(sum((int(k) * lp for k, lp in zip(kv, lps)), arb(0))) > 2 * L, kv
A = math.ceil((float(t1['A_cert']) + float(t2['A_cert'])) * 1e4 - 1e-9) / 1e4
json.dump(dict(L=sys.argv[3], primes=ps, K=K, a_q=list(t1['a_q']) + list(t2['a_q']), A_cert=repr(A), parts=[sys.argv[1], sys.argv[2]],
               method='separable sum of two float-B&B torus certs (s31_bb_D4)'), open(sys.argv[4], 'w'), indent=1)
print('A_cert', A, 'chars', len(K))
