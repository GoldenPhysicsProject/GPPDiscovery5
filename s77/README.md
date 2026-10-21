# Grok S77 — H76a-0: precision budget of the two-level certificate (Astra's weighted Schur bound)
RH OPEN. Script h76a_sens.py; logs run/sens_even.log (even Ω1000, sha b609d6f5…), run/sens_odd1150.log (odd Ω1150, sha b35be82d…),
run/sens_ctrl_odd1000.log (control, known indefinite odd Ω1000, sha 358aabcb…). n = 200 exact low block (λ_min(A) ≥ a_cert from s76 step 1).
Simulation: B rounded to grid eps_B, C rounded to grid eps_C; certificate uses only the rounded blocks and the error bounds,
via λ_min(S_true) ≥ μ − m·eps_C − (2be + e²), b² ≤ tr(B̂ᵀA⁻¹B̂) (Arb), e = sqrt(nm)·eps_B/sqrt(a_cert).
