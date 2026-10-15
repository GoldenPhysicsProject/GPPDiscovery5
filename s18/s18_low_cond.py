"""S18 claims 4 + 5: Ritz (NUMERICAL, outer/upper) lambda_min of the TRUE completed real-g Weil form at h=h*=log3/2 for every
residual primitive character class (f, parity, Re eta(2)) plus zeta and F_{1/4} (S4/S17 convention) controls.
 Q_eta = K_e + N (T_e(2h) + log f - c_e) - 2 (log2)/sqrt2 Re eta(2) T_2     (only n=2 active on (0,h*]; F(log3)=0)
Each Ritz value is accompanied by the Arb Rayleigh quotient of an explicit rational g (CERTIFIED upper bound on inf Q/N),
F_g(log2)/N of that minimiser, and the analytic CERTIFIED lower bound LB_C = E_comb(h*) + log f - c_e - (log2)/sqrt2 |Re eta(2)|.
usage: python s18_low_cond.py M1 M2 ..."""
import sys, json
from s18_ritz_common import *
from s18_common import HERE, forms
C0, C1 = consts(); EC = {0: E_comb(HSTAR, 0), 1: E_comb(HSTAR, 1)}; T2h = {0: Tarb(2*HSTAR, 0), 1: Tarb(2*HSTAR, 1)}
RES = {0: [5, 7, 8, 9, 11, 12, 13], 1: [3, 4, 5, 7]}
classes = {}
for e in (0, 1):
    for f in RES[e]:
        for l, par, o, e2, _ in primitive_chars(f):
            if par != e: continue
            re = e2.real; key = (f, e, round(float(re.mid()), 6))
            classes.setdefault(key, []).append(l)
# exact Re eta(2) for the classes (roots of unity): rational or cos(2 pi k/o)
EXACT = {-1.0: arb(-1), 0.0: arb(0), 1.0: arb(1), -0.5: arb(-1)/2, 0.5: arb(1)/2,
         0.309017: (2*arb.pi()/5).cos(), -0.809017: (4*arb.pi()/5).cos()}
rows = []
for M in map(int, sys.argv[1:]):
    B = build(M, HSTAR); N, T2 = B["N"], B["T"][2]
    F = forms(B, HSTAR)
    jobs = []
    for (f, e, r), labs in sorted(classes.items()):
        Ke, Ce = (B["K0"], C0) if e == 0 else (B["K1"], C1)
        reA = EXACT[r]
        Q = Ke + N*(T2h[e] + arb(f).log() - Ce) - T2*(2*LOG2S2*reA)
        lb = EC[e] + arb(f).log() - Ce - LOG2S2*abs(reA)
        jobs.append((f"f={f} par={e} Re eta(2)={r} labels={labs}", Q, lb))
    jobs.append(("zeta (control)", F["zeta"], None)); jobs.append(("F_1/4 (S4 conv, control)", F["F_theta"], None))
    for nm, Q, lb in jobs:
        for sn, I in SECT(M):
            lam, cert, fx = ritz_min(Q, N, I, extra=T2)
            row = dict(M=M, form=nm, sector=sn, ritz_min=mp.nstr(lam, 12), cert_upper=cert.str(12), F_log2_over_N=fx.str(8),
                       analytic_LB=(lb.str(10) if lb is not None else "n/a"))
            rows.append(row); print(json.dumps(row), flush=True)
json.dump(rows, open(os.path.join(HERE, "low_cond_" + "_".join(sys.argv[1:]) + ".json"), "w"), indent=1)
