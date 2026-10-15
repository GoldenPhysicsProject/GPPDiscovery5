"""S27: rationalize an LP multiplier at h* = (1/2) log 3 EXACTLY.
h is NOT rationalized: stored as h_spec = "log3/2". Only a[], y[] (den 1e15), D, c and the W-knot OFFSETS
w_off_q (z_j = 2h + w_off_j, w_off_0 = 0) are rational. RH OPEN.
usage: lp.json p/q out.json"""
import json, sys, numpy as np
from fractions import Fraction as F
d = json.load(open(sys.argv[1])); cp, cq = map(int, sys.argv[2].split("/")); den = 10**15
h = np.log(3)/2; assert abs(d["h"] - h) < 1e-15
D = F(d["D"]).limit_denominator(1000)
off = [F(x - 2*h).limit_denominator(10**6) for x in d["pl_knots"]]; off[0] = F(0)
assert all(abs(float(o) + 2*h - x) < 1e-12 for o, x in zip(off, d["pl_knots"]))
out = dict(h_spec="log3/2", D_q=[D.numerator, D.denominator], den=den, c=[cp, cq], Xi=d["Xi"],
           a_num=[round(x*den) for x in d["a"]], y_num=[round(x*den) for x in d["v_pl"]],
           w_off_q=[[q.numerator, q.denominator] for q in off], source=sys.argv[1], c_LP=d["c_LP"])
json.dump(out, open(sys.argv[3], "w"), indent=1); print("wrote", sys.argv[3])
