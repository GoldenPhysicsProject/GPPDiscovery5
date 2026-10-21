"""Grok S67-resumed: NUMERICAL lambda_min of a dumped midpoint matrix T (inverse power iteration on the shifted
midpoint matrix, flint arb_mat point arithmetic). Used only to pick l0 for the PD certificates; not a certificate.
usage: lam_from_dump.py dump.json [prec]"""
import sys, json
from flint import arb, arb_mat, ctx
ctx.prec = int(sys.argv[2]) if len(sys.argv) > 2 else 512
S = json.load(open(sys.argv[1])); n = len(S)
A = arb_mat([[arb(x) for x in r] for r in S])
Ai = A.solve(arb_mat([[arb(1) if i == j else arb(0) for j in range(n)] for i in range(n)]), algorithm='approx').mid()
v = arb_mat([[arb(1) / (i + 1)] for i in range(n)])
for it in range(60):
    u = (Ai * v).mid(); nu = sum((u[i, 0] ** 2 for i in range(n)), arb(0)).sqrt()
    lam = sum((v[i, 0] ** 2 for i in range(n)), arb(0)) / sum((u[i, 0] * v[i, 0] for i in range(n)), arb(0))
    v = arb_mat([[(u[i, 0] / nu).mid()] for i in range(n)])
print(json.dumps(dict(dump=sys.argv[1], K=n, lam_min_est=lam.mid().str(12, radius=False))))
