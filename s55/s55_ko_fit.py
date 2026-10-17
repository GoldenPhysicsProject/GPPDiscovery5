"""Grok S55 H55b: knockout thresholds L_p* (first L where W_zeta minus the prime-p term has a NEGATIVE Ritz value,
M=56 S14/S44 exact basis) -> delta_p = L_p* - (1/2)log p; fit power laws in p and in w_p = log p / sqrt p.
Intervals: (last positive h, first negative h]; Ritz >= inf so the right end is a rigorous-up-to-numerics UPPER bound.
NUMERICAL data, fitted law = CONJECTURE/HEURISTIC.  usage: s55_ko_fit.py  (reads ../s54/h54b and h55b)"""
import json, glob, math, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
rows = []
for f in glob.glob(os.path.join(here, '..', 's54', 'h54b', 'ko_M56*.jsonl')) + glob.glob(os.path.join(here, 'h55b', 'ko_*.jsonl')):
    for l in open(f):
        try: r = json.loads(l)
        except Exception: continue
        if r.get('M') != 56: continue
        rows.append(r)
res = {}
for p in (2, 3, 5, 7, 11, 13, 17, 19):
    hp = 0.5 * math.log(p); res[p] = {}
    for s in ('even', 'odd'):
        pos, neg = [], []
        for r in rows:
            v = r[s].get(f'no{p}')
            if v is None: continue
            h = float(r['h'])
            if h <= hp: continue
            (neg if float(v[0]) < 0 else pos).append(h)
        if not neg: continue
        hi = min(neg); lo = max([h for h in pos if h < hi], default=hp)
        res[p][s] = dict(lo=lo, hi=hi, d_lo=lo - hp, d_hi=hi - hp, d_mid=0.5 * (lo + hi) - hp)
out = {'label': 'NUMERICAL (Ritz M=56); fits CONJECTURE', 'half_log_p': {p: 0.5 * math.log(p) for p in res}, 'delta': res, 'fits': {}}
for s in ('odd', 'even'):
    P = [p for p in res if s in res[p]]
    d = np.array([res[p][s]['d_mid'] for p in P]); w = np.array([math.log(p) / math.sqrt(p) for p in P]); pp = np.array(P, float)
    for nm, x in (('p', pp), ('w_p', w), ('log_p', np.log(pp))):
        b, a = np.polyfit(np.log(x), np.log(d), 1); pred = np.exp(a) * x ** b
        out['fits'][f'{s}_vs_{nm}'] = dict(exponent=float(b), coef=float(math.exp(a)), rms_log_resid=float(np.sqrt(np.mean((np.log(d) - np.log(pred)) ** 2))), primes=P)
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(here, 'h55b', 'ko_fit.json'), 'w'), indent=1)
