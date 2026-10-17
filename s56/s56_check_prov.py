"""Grok S56: validate the provisional-lift assembly. PASS iff the B&B-certified merged lift has identical primes/K/a_q
to the provisional JSON used in the assembly and A_cert(certified) <= A_cert(provisional) (a larger A in the assembly only
lowers beta and the lift constant, so the assembly with A_prov is a valid lower bound). RH OPEN.
usage: python s56_check_prov.py certified.json provisional.json"""
import json, sys
c = json.load(open(sys.argv[1])); p = json.load(open(sys.argv[2]))
assert c['primes'] == p['primes'], 'primes differ'; assert c['K'] == p['K'], 'K differs'; assert c['a_q'] == p['a_q'], 'a_q differ'
ac, ap = float(c['A_cert']), float(p['A_cert'])
print('certified A_cert', ac, 'provisional', ap, 'PASS' if ac <= ap else 'FAIL (assembly invalid, rerun with A_cert)')
