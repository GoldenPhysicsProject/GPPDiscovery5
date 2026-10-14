"""CERTIFIED (Arb ball arithmetic, python-flint) Li coefficients of F_theta(s) = xi(s-theta) xi(s+theta).

lambda_n = n [z^n] log F_theta(1/(1-z)). We build the exact formal power series u = z/(1-z) and
evaluate log xi(1 -/+ theta + u) with acb_series (log, lgamma, zeta, all rigorous ball arithmetic).
xi(s) = s(s-1)/2 pi^{-s/2} Gamma(s/2) zeta(s). Every coefficient is a rigorous enclosure; a sign is
certified when the ball excludes 0.
Run: python experiments/li_certified/li_arb.py L prec [theta_num theta_den]
"""
import os, sys, time, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from flint import acb_series, acb, arb, ctx, fmpq
from common import write_result
L = int(sys.argv[1]); ctx.prec = int(sys.argv[2]); ctx.cap = L + 1
tn, td = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (1, 4)
th = arb(fmpq(tn, td))
t0 = time.time()
u = acb_series([0] + [1] * L)
lpi = acb(arb.pi()).log()
def logxi(x):
    return (x * (x - 1) / 2).log() - (x / 2) * lpi + (x / 2).lgamma() + x.zeta().log()
c = (logxi(acb(1 - th) + u) + logxi(acb(1 + th) + u)).coeffs()
lam = {n: n * c[n] for n in range(1, len(c))}
neg = [n for n, v in lam.items() if v.real < 0]          # certified negative (ball entirely < 0)
undecided = [n for n, v in lam.items() if not (v.real < 0) and not (v.real > 0)]
first_neg = min(neg) if neg else None
print(f"theta={tn}/{td} L={L} prec={ctx.prec} time={time.time()-t0:.1f}s")
print("certified negative count:", len(neg), "first:", first_neg, "undecided:", undecided[:20])
for n in sorted(set([1, 100, 1000] + ([first_neg - 1, first_neg] if first_neg else []) + [L])):
    if n in lam: print(n, lam[n].real.str(15, radius=True))
tag = f"{tn}_{td}_L{L}"
write_result(os.path.join(os.path.dirname(__file__), f"result_{tag}.json"),
  {"experiment": "li_arb", "label": "CERTIFIED sign (Arb balls)", "theta": f"{tn}/{td}", "L": L,
   "first_certified_negative_n": first_neg, "all_positive_below": (first_neg or L+1) - 1,
   "undecided": undecided,
   "lambda_balls": {str(n): lam[n].real.str(20, radius=True) for n in lam if n % 50 == 0 or n in (1, first_neg, (first_neg or 2)-1)}},
  kind="certified", method="python-flint acb_series ball arithmetic", dps=int(ctx.prec*0.30103))
