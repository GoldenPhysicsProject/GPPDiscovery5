"""S24: rationalise an S22 LP symbol cert + design a continuous PL extension (graded knots, honest tail) -> Arb-ready JSON."""
import json, sys, numpy as np
from fractions import Fraction as Fr
import s24_design as S
src, hs, H, N1, N2, g1, g2, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), sys.argv[8]
cfrac = float(sys.argv[9]) if len(sys.argv) > 9 else 1.0
d = json.load(open(src)); hq = Fr(hs); d["h"] = str(float(hq))
d["h_q"] = [hq.numerator, hq.denominator]; d["D_q"] = [1, int(round(1/d["D"]))]; d["Dd_q"] = [1, int(round(1/d["Dd"]))]
cq = Fr(d["c_target"]*cfrac).limit_denominator(10000)
if float(cq) > d["c_target"]*cfrac: cq = Fr(cq.numerator-1, cq.denominator)
d["c"] = [cq.numerator, cq.denominator]
zf = S.knots2(float(hq), H, N1, N2, g1, g2); zden = 10**9
zq = [Fr(round(z*zden), zden) for z in zf]; zq[0] = hq; zq[-1] = Fr(H).limit_denominator(1000)
zs = np.array([float(z) for z in zq])
o = S.design_tail(d, zs)
aden = 10**12; al = [int(round(x*aden)) for x in o["alpha"]]
y0 = np.sinh(float(hq)/2); y0n = int(round(y0*aden))
N = len(zs)-2
d["extension"] = dict(type="piecewise_linear_odd_continuous", H_q=[zq[-1].numerator, zq[-1].denominator], N=N+1,
    alpha_num=[y0n]+al, alpha_den=aden, z_num=[int(z*zden) for z in zq], z_den=zden,
    knots=dict(N1=N1, N2=N2, g1=g1, g2=g2), two_beta_design_float=o["two_beta"], mX=o["mX"], Xi_design=o["Xi"])
# N in verifier = number of free values incl. y at h; zs has N+2 entries -> here N+1 values + final 0 => len(z)=N+2 with N=len(alpha)
d["extension"]["N"] = len(d["extension"]["alpha_num"]) - 1
assert len(d["extension"]["z_num"]) == d["extension"]["N"]+2, (len(d["extension"]["z_num"]), d["extension"]["N"])
json.dump(d, open(out, "w")); print(out, "c=", cq, float(cq), "2beta_design=", o["two_beta"], "N=", d["extension"]["N"])
