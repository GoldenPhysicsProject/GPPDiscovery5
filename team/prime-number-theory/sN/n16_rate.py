# N.16.3 derived margin formula (CONJECTURE), evaluated against the N.15 Ritz upper bounds. No Ritz scans here.
# D(t*) = L t*/pi - Nsm(t*) - 1 at t* = 2 pi e^{2L} (per parity), Nsm = RvM smoothed count; c* = L t*;
# margin^sigma ~ exp(-pi^2 (2D - sigma)/log c*), sigma = 0 even, 1 odd; capped at O(1) when 2D - sigma <= 0.
import math, json
g = json.load(open(__file__.replace("n16_rate.py", "n16_zeros400.json")))
Nsm = lambda T: T/(2*math.pi)*math.log(T/(2*math.pi*math.e)) + 7/8
N15 = {"log2": (0.073, 0.5493), "log2+0.1": (0.0168, 0.2452), "log2+0.3": (4.5e-4, 0.0142), "log3": (2.5e-4, 1.34e-3)}
rows = [("log2", math.log(2)), ("log2+0.1", math.log(2)+0.1), ("log2+0.3", math.log(2)+0.3), ("log3", math.log(3)),
        ("log3.5", math.log(3.5)), ("log4", math.log(4)), ("log6", math.log(6)), ("log10", math.log(10))]
print(" 2L        mu     t*      N(t*)  Nsm    D(t*)   c*      even~      odd~       N.15 Ritz<= (even, odd)")
for name, twoL in rows:
    L = twoL/2; mu = math.exp(twoL); ts = 2*math.pi*mu; c = L*ts
    D = L*ts/math.pi - Nsm(ts) - 1
    f = lambda s: 1.0 if 2*D - s <= 0 else math.exp(-math.pi**2*(2*D - s)/math.log(c))
    Nt = sum(1 for x in g if x <= ts)
    print(f" {name:9s} {mu:5.2f} {ts:7.3f}  {Nt:4d}  {Nsm(ts):6.3f} {D:6.3f}  {c:6.2f}  {f(0):9.2e}  {f(1):9.2e}  {N15.get(name, ('-', '-'))}")
print("odd onset 2D = 1: mu =", 15/8 + 0.5, " 2L =", round(math.log(15/8 + 0.5), 4), " ell2 = 2L - log2 =", round(math.log(15/8 + 0.5) - math.log(2), 4))
print("even onset 2D = 0: mu =", 15/8, " 2L =", round(math.log(15/8), 4), "(< log 2: even sector is in the plunge from the p = 2 entry)")
