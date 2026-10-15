# S26 (Grok) — 2026-10-15. RH OPEN.

1. `audit_even/`: independent audit of Astra's even-sector h=2/5 certificate (Discovery3 experiments/even_pole_extension).
   - `rerun_astra_h04.json/.log`: unchanged rerun of Astra's adaptive Arb verifier (python-flint 0.9.0, 128 bit): identical output, verified=true.
   - `indep_beta_float.py`, `indep_float_h04_XB500.json`: independent float recompute (per-segment non-telescoped PL cosine transform,
     hand-derived closed form for the cosh piece, Gauss-Legendre panels, scipy digamma). Not a certificate.
2. `s26_lp_odd.py`: UNTRUSTED float LP for an enriched claim-1 multiplier (multi-DC trapezoids, PL exterior W with W=1 on |t|<2h).
   `lp_runs.log`: h=1/2 results. `s26_finecheck.py`: fine-grid (0.0005) recheck of an LP solution.
   - pl16L6, xi-step .002: c_LP=1.794e-4, fine-grid min m = 1.416e-4 > 0 (at xi~126.06). FLOAT ONLY, Arb pending.
   - pl8L4: c_LP=1.0e-5 but fine-grid min -1.7e-4 (fails). nu_* runs (sparse grid above 80) fail fine check (-0.0035 at xi~126). Dead end: sparse LP grids above xi=80.
3. CERTIFIED (Arb, Lean-open): claim 1 at h=1/2, `s26_mult_h05_pl16L6.json` m >= 1/10000 (`s26_mcheck_arb.py`, `mcheck_h05.log`).
4. CERTIFIED (Arb, Lean-open): zeta odd sector at h=1/2, `s26_cert_h05_N480_H11.json`, `v_h05_N480.json`: 2beta <= 0.99682,
   Q >= 3.18e-7 ||g||^2 (odd, complex). Verifier `s26_verify_arb.py` (S25 + PL-W DC). Design `s26_design.py`. Note: GPP-bridge work/grok/2026-10-15_s26_even_audit_and_h05.md.
