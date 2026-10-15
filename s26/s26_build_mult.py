"""Rationalize S26 LP multiplier (PL exterior) to denominator 1e15. c target given as p/q."""
import json, sys
from fractions import Fraction as F
d = json.load(open(sys.argv[1])); cp, cq = map(int, sys.argv[2].split("/")); den = 10**15
h = F(d["h"]).limit_denominator(1000); D = F(d["D"]).limit_denominator(1000)
z = [F(x).limit_denominator(10**6) for x in d["pl_knots"]]
assert all(abs(float(a)-b) < 1e-14 for a, b in zip(z, d["pl_knots"]))
out = dict(h_q=[h.numerator, h.denominator], D_q=[D.numerator, D.denominator], den=den, c=[cp, cq], Xi=d["Xi"],
           a_num=[round(x*den) for x in d["a"]], y_num=[round(x*den) for x in d["v_pl"]],
           z_q=[[q.numerator, q.denominator] for q in z], source=sys.argv[1], c_LP=d["c_LP"])
json.dump(out, open(sys.argv[3], "w"), indent=1); print("wrote", sys.argv[3])
