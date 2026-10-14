"""Golden ratio in the Davenport-Heilbronn coefficient. Exact identities + CERTIFIED Arb check.
chi mod 5 with chi(2)=i. tau(chi) = -2 sin36 + 2i sin72 (exact), eps(chi)=tau/(i sqrt5) = (phi+i)/sqrt(phi+2).
DH functional equation forces (1-i kappa)/(1+i kappa) = eps(chibar), so kappa = tan(arctan(1/phi)/2)
= sqrt(phi+2) - phi = 5^{1/4} sqrt(phi) - phi = (sqrt(10-2sqrt5)-2)/(sqrt5-1).
Also class-number formula for Q(sqrt5): L(1,chi_5) = 2 log(phi)/sqrt5.
phi here is FORCED BY THE CHOICE q=5 (Q(zeta_5) contains Q(sqrt5)); not forced by RH/fakes in general.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from flint import arb, acb, ctx
from common import write_result
ctx.prec = 400
s5 = arb(5).sqrt(); phi = (1 + s5) / 2
z = lambda k: acb.exp_pi_i(acb(arb(2 * k) / 5))
chi = {1: acb(1), 2: acb(0, 1), 4: acb(-1), 3: acb(0, -1)}
tau = sum((chi[n] * z(n) for n in chi), acb(0))
eps = tau / (acb(0, 1) * s5)
eps_phi = acb(phi, 1) / (phi + 2).sqrt()
kappa_def = ((10 - 2 * s5).sqrt() - 2) / (s5 - 1)
kappa_phi = (phi + 2).sqrt() - phi
k2 = arb(5).root(4) * phi.sqrt() - phi
a = acb(1, -kappa_phi)
fe = a / a.conjugate() - eps.conjugate()
checks = {"eps - (phi+i)/sqrt(phi+2)": (eps - eps_phi), "kappa_def - (sqrt(phi+2)-phi)": acb(kappa_def - kappa_phi),
          "kappa - (5^{1/4}sqrt(phi)-phi)": acb(kappa_def - k2), "(1-ik)/(1+ik) - eps(chibar)": fe}
for k, v in checks.items(): print(k, v.str(5, radius=True), "contains0:", v.contains(0))
print("kappa =", kappa_phi.str(30))
write_result(os.path.join(os.path.dirname(__file__), "result.json"),
  {"experiment": "dh_phi", "label": "PROVED LEMMA (elementary; classical DH constant rewritten) + CERTIFIED ball check",
   "kappa": kappa_phi.str(40), "checks": {k: v.str(5, radius=True) for k, v in checks.items()},
   "phi_status": "comes from the choice q=5"}, kind="certified", method="python-flint arb, 400 bits")
