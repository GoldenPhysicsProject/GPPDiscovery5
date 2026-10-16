"""Grok S35: structural audit of a merged (separable) torus certificate: every invisible character must be supported
in exactly one prime group, and K/a_q must be exactly the concatenation of the parts (so sup = sum of part sups).
Usage: python s35_check_merge.py merged.json part1.json[:offset] part2.json[:offset] ...   RH OPEN."""
import sys, json
m = json.load(open(sys.argv[1])); n = len(m['primes']); K = []; A = []; groups = []
for spec in sys.argv[2:]:
    fn, off = (spec.split(':') + ['0'])[:2]; off = int(off); p = json.load(open(fn))
    k_p = len(p['primes']); assert m['primes'][off:off + k_p] == p['primes'], (fn, p['primes'])
    groups.append(set(range(off, off + k_p)))
    for kv, a in zip(p['K'], p['a_q']):
        v = [0] * n; v[off:off + k_p] = kv; K.append(v); A.append(a)
assert m['K'] == K and m['a_q'] == A, 'merged K/a_q differ from concatenation of parts'
for kv in m['K']:
    sup = {i for i, k in enumerate(kv) if k}
    assert sum(1 for g in groups if sup & g) == 1, ('character straddles groups', kv)
print(sys.argv[1], 'OK: concatenation of parts, all', len(K), 'characters group-local; A_cert =', m['A_cert'])
