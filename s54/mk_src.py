"""Grok S54: turn a bundle-LP JSON into the B&B source JSON (adds per_prime = sum of visible coefficients, lam_max=None).
usage: mk_src.py bundle.json out_src.json"""
import sys, json
d = json.load(open(sys.argv[1])); k = list(d)[0]; v = d[k]
v['per_prime'] = sum(c for (_, _, c) in v['vis']); v['lam_max'] = None
json.dump({k: v}, open(sys.argv[2], 'w')); print(k, 'A_true', v['A_true'], 'chars', len(v['K']))
