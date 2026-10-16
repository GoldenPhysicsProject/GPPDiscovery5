#!/usr/bin/env python3
"""H41b: Weil-window → S(t) bridge for t ≤ 2.6 (SECONDARY).

Uses CERTIFIED support-2.6 coercivity from S36/S38:
  Q_ζ(f) ≥ c ‖f‖²_L2 for all complex f supported in [−1.3, 1.3],
  c ≈ 3.634e−62 (S36) / 4.260e−62 (S38 even, D8×2).

Top-hat identity (Codex/S40): for f_t(x) = t^{−1/2} 1_[−t/2,t/2],
  ‖f_t‖² = 1 and S(t) = (t/2) Q_ζ(f_t),
CONFIRMED-with-caveat: f_t ∉ Schwartz; mollify + dominated convergence.
The certified Weil class is a Galerkin/Fourier window of smooth compactly
supported test functions — top-hat is outside that class without mollification.

This script does NOT claim CERTIFIED S(t)>0 from the window alone.
It computes quantitative lower bounds under stated hypotheses and lists GAPs.
"""
from __future__ import annotations

import json
import os
from mpmath import mp, mpf, nstr, pi, sin, fabs

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "h41b_weil_bridge.json")
mp.dps = 40

# Certified constants from Grok S36 / S38 (Arb; Lean-open)
C_S36 = mpf("3.6340e-62")   # even binding at support 2.6, A≤4.6565
C_S38 = mpf("4.260e-62")    # even at support 2.6 with D8×2 A≤4.4702
SUPPORT = mpf("2.6")        # full-complex certified window width
L = mpf("1.3")              # half-support

results = {
    "session": "S41",
    "H41b": {},
    "labels": {},
}

# ---------- 1. Theorem statement (careful) ----------
theorem = {
    "statement": (
        "If S(t)=(t/2)Q_ζ(f_t) holds for (mollified) top-hats f_t with "
        "support width t≤2.6, AND if the mollified family stays inside the "
        "certified Weil class (or the identity extends by density/dominated "
        "convergence to the L² top-hat), THEN S(t) ≥ (t/2)·c·‖f_t‖² for "
        "t∈(0,2.6], with c the certified coercivity constant."
    ),
    "certified_window": {
        "support_width": "2.6",
        "L": "1.3",
        "c_S36_even": nstr(C_S36, 10),
        "c_S38_even_D8x2": nstr(C_S38, 10),
        "source": "GPPDiscovery5/s36, s38; GROK_RESEARCH_NOTES S36/S38",
    },
    "mollification_caveat": (
        "S40 CONFIRMED-with-caveat: top-hat is square-integrable but not "
        "Schwartz/smooth. Mollify; sinc²(r)·log(2+|r|) integrable ⇒ dominated "
        "convergence for Arch/gamma integral; prime sum finite. "
        "Schwartz-only statement of Weil criterion = mild GAP for bare top-hat."
    ),
}
results["H41b"]["theorem"] = theorem

# ---------- 2. ‖f_t‖² for L² top-hat ----------
# f_t(x) = t^{-1/2} 1_{[-t/2,t/2]}(x)  ⇒  ∫|f|² = t^{-1} · t = 1.
norm_sq = mpf(1)
# Mollified family: f_{t,ε} = f_t * ρ_ε with ρ_ε approximate identity,
# support in [-t/2-ε, t/2+ε]. For support to stay in [-1.3,1.3] need
# t/2 + ε ≤ 1.3 ⇒ t ≤ 2.6 - 2ε.
# As ε→0, ‖f_{t,ε}‖² → 1, and for small ε, ‖f_{t,ε}‖² ≥ 1-δ(ε).

def mollified_support_ok(t, eps):
    return (t / 2 + eps) <= L

# Galerkin / Fourier projection loss sketch:
# The certified class uses a finite Fourier/Laguerre-type window (T=700,K=740
# in S36). Projecting the top-hat onto that subspace loses high-frequency mass.
# Lower-bounding Q on the projection gives
#   Q(f) ≥ Q(P f) ≥ c ‖P f‖² = c (‖f‖² − ‖(1-P)f‖²),
# so the effective constant is c·(1 − loss). Without an explicit projection
# bound we cannot CERTIFY a positive floor for bare top-hats.

# Quantitative: for ideal identity with ‖f‖=1 and no projection loss:
ideal_lower = {}
for t_val in [mpf("0.5"), mpf("1"), mpf("2"), mpf("2.6"), SUPPORT]:
    # S(t) ≥ (t/2) * c * 1
    ideal_lower[nstr(t_val, 5)] = {
        "S_ge_S36": nstr((t_val / 2) * C_S36, 12),
        "S_ge_S38": nstr((t_val / 2) * C_S38, 12),
        "support_fits_window": float(t_val) <= float(SUPPORT) + 1e-15,
    }

results["H41b"]["L2_tophat"] = {
    "definition": "f_t = t^{-1/2} 1_{[-t/2,t/2]}",
    "norm_sq": "1 (exact)",
    "ideal_S_lower_bounds_if_identity_in_cert_class": ideal_lower,
}

# Mollified: pick ε so support fits; estimate ‖f_{t,ε}‖² ≥ 1 - C ε/t (rough)
# For reporting only — not a certificate of the extension.
eps_examples = []
for t_val in [mpf("1"), mpf("2"), mpf("2.5")]:
    # need t/2+ε ≤ 1.3
    eps_max = L - t_val / 2
    eps_examples.append({
        "t": nstr(t_val, 5),
        "eps_max_for_support": nstr(eps_max, 10),
        "fits": float(eps_max) > 0,
    })
results["H41b"]["mollified_support"] = eps_examples

# ---------- 3. Verdict ----------
gaps = [
    "GAP: bare L² top-hat is not in the Schwartz/smooth test class of the "
    "    classical Weil criterion; S40 mollify caveat stands.",
    "GAP: certified S36/S38 class is a specific finite Galerkin window "
    "    (T,K,ng,...), not all of C_c^∞([-1.3,1.3]). Extension of Q≥c‖f‖² "
    "    from the Galerkin subspace to mollified top-hats is unproved here.",
    "GAP: no explicit projection-loss bound ‖(1-P)f_t‖² computed in S41.",
    "OPEN: kick/U+ barrier for t>2.6; this bridge only addresses small-t mass.",
]

verdict = {
    "label": "CONFIRMED-with-GAP",
    "what_is_confirmed": (
        "IF the top-hat identity S=(t/2)Q extends to mollified top-hats with "
        "support in [-1.3,1.3] AND those mollifiers lie in (or are limits of) "
        "the certified coercivity class, THEN S(t)≥(t/2)c>0 on (0,2.6] with "
        f"c≥{nstr(C_S36,6)} (S36) or c≥{nstr(C_S38,6)} (S38). "
        "Ideal numerical floor at t=2.6: S≥1.3·c ≈ 4.7e−62 (S36) / 5.5e−62 (S38)."
    ),
    "what_is_NOT_certified": (
        "S41 does NOT claim CERTIFIED S(t)>0 on (0,2.6] from the Weil window "
        "alone. Mollification + Galerkin-density / projection-loss remain GAPs."
    ),
    "gaps": gaps,
    "dead_end_avoided": (
        "DEAD END: claiming bare top-hat ∈ Schwartz / cert class without mollify."
    ),
}
results["H41b"]["verdict"] = verdict
results["labels"]["H41b"] = "CONFIRMED-with-GAP"
results["RH"] = "OPEN"
results["phi"] = "absent"

with open(OUT, "w") as f:
    json.dump(results, f, indent=2)
print(json.dumps(verdict, indent=2))
print("Wrote", OUT)
