# S22 — ζ odd sector on the first prime window: pole Schur step + odd-only DC lift (Arb certificate at h = 2/5)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH OPEN.**

- `s22_odd_susc.py M h...` (NUMERICAL, Ritz/Galerkin, S17 builder): λ_min of Q_ζ and of the pole-free form B on odd g, and the
  Galerkin susceptibility χ = 2vᵀB⁻¹v, v = sinh(x/2) (GPPVerify `OddCollisionTrace.galerkin_criterion`). Output `out_susc_M24.json`.
- `s22_lp_odd.py` (UNTRUSTED float LP design): positive symbol m = σ_A + C·w_DC − Σ a_k cos(u_k ξ), u_k = 2h + k/20, where w_DC is the
  Fourier transform of a trapezoid equal to 1 on |t| ≤ 2h. Both added families are invisible to odd g on (−h,h). `save h` writes
  `s22_cert_h<h>.json` (rational a_k, C, target c).
- `s22_verify_arb.py` (RIGOROUS, Arb 128-bit): Claim 1 m ≥ c on ℝ (S20 cell method + DC derivative bounds + digamma-monotone tail);
  Claim 2 an upper bound for β = (1/2π)∫|v̂|²/m with v = sinh(x/2)·1_(−h,h). Output `s22_verify_results.json`.
- `s22_symbol_check.py` (NUMERICAL): S17 matrix (Q_ζ − P)_jk equals (1/2π)∫ĝ_j ĝ_k σ_A and (1/2π)∫ĝ_j ĝ_k m to 9 digits
  (odd basis, h = 2/5, ξ-grid 0.002 up to 800). `out_symbol_check.txt`.
- `s22_beta_float.py`: float preview of β.

Result at h = 2/5: m ≥ c = 44/3917 verified; 2β ≤ 0.82970 < 1 verified; hence Q_ζ(g) ≥ (1−2β)c‖g‖² ≥ 0.001913‖g‖² for odd g
(real or complex) in C_c((−2/5, 2/5)). h = 9/20: claim 1 passes (c = 38/4235), but the un-extended v gives 2β ≈ 1.50 (fails);
it needs an optimized extension of v outside the window. See the bridge note for the proof.
