# Grok S76 — H75b two-level (low-block exact LDLᵀ + Schur complement) certificate, tested on the support-2.84 dumps
RH OPEN. Dumps: ../s71/run/T_L7150_m12_O1000_p0.json (even, sha b609d6f5…), T_L7150_m12_O1150_p1.json (odd, b35be82d…),
control T_L7150_m12_O1000_p1.json (odd, known indefinite, sha 358aabcb…).
- h75b_two_level.py: sweep n ∈ {1180,1150,1100,1000,800,400,200}: certified g_C (Gershgorin of C) and ‖B‖_F, numerical λ_min(A),
  Schur complement S diagnostics, scaled variant. Log sweep_even_O1000.log / sweep_even_O1000.json.
- h75b_cert.py: rigorous Schur S (Arb solve) + Gershgorin on S (n = 1180), run/cert_even_n1180.log.
- h75b_cert2.py: rigorous Schur S + float Cholesky of S − tI with exact Arb residual (n = 200), run/cert2_n200.log, run/cert2_more.log.
- subdump.py: leading n×n block for ../s71/s37_ldlt_exact.py (step 1), logs run/ldlt_A_*.log.
