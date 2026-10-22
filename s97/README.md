# Grok S97 (2026-10-22): H96a, the discrete Kreĭn/de Branges chain of the nested window family. RH OPEN.

Status: **HIGH-PRECISION NUMERICAL** (mpmath dps 60). It is not interval arithmetic. Every PD/fail is cross-checked by mp Cholesky of the full Toeplitz matrix, both before and at the failure index.

## Model (`chain.py`)
- The windows [0, T] are nested. The basis is step functions χ_k on cells [kh, (k+1)h].
- The Weil form is translation invariant, so the Gram matrix is Toeplitz: c_m = W(even part of tent_h(· − mh)).
- Levinson–Durbin gives the Schur/Verblunsky coefficients κ_k and pivots P_k = P_{k−1}(1 − κ_k²). T_n ≻ 0 ⇔ c_0 > 0 and |κ_k| < 1 for all k < n.
- κ_k/h → Kreĭn's continuous reflection function as h → 0. That function determines the de Branges Hamiltonian H(t), so this is the discrete canonical system for the window family.
- **Calibration (CHECKED):** with n = 2 cells of width L, the model equals exactly 2 × the s87/hp_forms N = 1 values Q_e and Q_o for ζ, DH and F_0.1 (15 digits, L = 0.5). The factor 2 is the overall normalisation only.

## Runs
- `run.py` (h 0.1, 0.02) and `run2.py` (h 0.02, 0.01) cover T ≤ 3 (our L ≤ 1.5).
- Forms: ζ, DH, F_0.1, F_0.05. Λ-controls: drop p = 2, 3, 5, 7. Pole control: pole weight × 0.9.
- `summarize.py` writes `summary.txt`.

## Results
1. **Sign content (PROVED, standard: Szegő/Schur/Kreĭn).** "H(t) ⪰ 0 on [0, T]" ⇔ |κ_k| < 1 for all k ⇔ every nested window form is PD. So the chain carries exactly window positivity and no extra sign.
2. **Locality of the chain (CHECKED from structure; seen NUMERICALLY).**
   - Each c_m receives prime atoms only where |mh − log n| < h, so κ jumps only at lags ≈ log n.
   - ζ's top |Δ²κ| lags (h 0.01) are 1.08, 1.59, 2.18, 2.38, 2.54, 2.81 and 2.92, which match log 3, 5, 9, 11, 13, 17 and 19.
   - DH's are 1.38, 1.77, 2.19, 2.62 and 2.75, which match log 4, 6, 9, 14 and 16 (its own Dirichlet support).
   - This is automatic, and it is Suzuki's Σ_{n ≤ e^t} structure.
3. **Λ-defects fail locally (NUMERICAL, converged in h).**
   - Dropping p fails at T_fail − log p = 0.057, 0.041, 0.031 and 0.024 for p = 2, 3, 5, 7 (h 0.01; h 0.02 gives 0.067, 0.041, 0.031, 0.034).
   - Failure cannot come before log p, since the forms are identical below it. That it comes so soon after is the S90 Λ-rigidity again.
4. **Analytic fakes fail non-locally (NUMERICAL).**
   - F_0.1 is PD up to T = 3 at h 0.02, and fails at T = 2.03 at h 0.01. That point lies between log 7 and log 8 and is not converged in h.
   - The Galerkin onset (S96) is T = 2L ≈ 1.0, between log 2 and log 3. H₋₄'s Galerkin onset T ∈ (0.9, 1.0] contains no log n either. H₋₂'s bracket (1.3, 1.4] contains log 4 = 1.386, so it is inconclusive on the 0.1 grid.
   - So "indefiniteness of H(t) is confined to prime slabs" is FALSE. Analytic deformations fail globally.
5. **DH is PD throughout T ≤ 3** (h 0.01, Cholesky), as in Galerkin. The chain does not rule it out. Its κ signs differ from ζ's from T 0.67 on, but that reflects its different arch kernel and coefficients, not a sign separator.
6. **Resolution caveat.** Step functions cannot reach the e^{−cμ} near-null directions (S91). The step chain sees F_0.1 only at h ≤ 0.01, and cannot see H₋₂ at all. A piecewise-polynomial block chain would be needed for those, but item 1 makes it pointless as a separator.

**Verdict: DEAD END** as an RH separator. The chain only restates nested positivity, and its locality is automatic. Prior art:
- Suzuki 1204.1827 (canonical system for Θ_ω; RH ⇔ PSD Hamiltonian family, constructed for ω > 1, via Burnol's method);
- Suzuki 2206.03682 and 2209.04658 (screw function and line; t-local Σ_{n ≤ e^t} Λ(n)/√n terms);
- Lagarias (de Branges spaces and ζ, 2006);
- Burnol (Sonine spaces);
- Kreĭn's continuous analogue of orthogonal polynomials.

What it rules out: Λ-defects (local, just after log p), the pole-weight control and F_0.1 (non-locally, at fine h). It does not rule out DH, and H₋₂ is not resolved by steps.
