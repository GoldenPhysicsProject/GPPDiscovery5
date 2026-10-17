"""Grok S56: build the PROVISIONAL merged torus JSON for the buffered L=69/50 lift (same K and a_q = '%.9e' rounding as
s34_bb_D4 + s34_merge_torus2), with A_cert set to a provisional ceiling. The assembly result is valid iff the B&B-certified
merged A_cert <= this ceiling AND K/a_q are identical (checked by s56_check_prov.py after the B&B finishes). RH OPEN."""
import json, sys
A_PROV = sys.argv[1] if len(sys.argv) > 1 else '4.9610'
p1 = json.load(open('lift/bundle_L138_D9_p235_buf.json'))['L1.38_D9']; p2 = json.load(open('lift/bundle_L138_D8_p71113_buf.json'))['L1.38_D8']
n1, n2 = len(p1['primes']), len(p2['primes'])
K = [list(k) + [0] * n2 for k in p1['K']] + [[0] * n1 + list(k) for k in p2['K']]
a_q = ['%.9e' % x for x in p1['a_lp']] + ['%.9e' % x for x in p2['a_lp']]
json.dump(dict(L='69/50', primes=p1['primes'] + p2['primes'], K=K, a_q=a_q, A_cert=A_PROV, parts=['PROVISIONAL'],
               method='PROVISIONAL: A_cert pending s34_bb_D4/s35_interval_bb B&B of both parts (s56/lift)'),
          open('run138/torus_L1.38_buf_PROVISIONAL.json', 'w'), indent=1)
print('chars', len(K), 'A_prov', A_PROV)
