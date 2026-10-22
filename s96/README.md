# Grok S96 (2026-10-22): H95a, ‖M‖ ≤ 1 as an L → 2L doubling step. RH OPEN.

Status: **HIGH-PRECISION NUMERICAL** (mpmath dps 90 forms, Legendre N 24 on the half-line s ∈ [0, L], L 0.3–1.4). Every sign is from mp Cholesky (PD) or LDLᵀ (n₋, matching eigen counts). Failing controls are F_0.1, H₋₂ and H₋₄. Not interval arithmetic.

Run: `python h95a.py` (16 L values, 6 processes, ~1.5 min) → `h95a.json`, `h95a.log`. Then `python summarize.py` → `summary.txt`. Needs `forms_cache.pkl` → `../s91/forms_cache.pkl` (a symlink, not committed; DH/F_θ at the missing L are built on the fly, ~8 s each).

## Setup
A = (Q_e + Q_o)/4 is the form on functions supported in [0, L]. B = (Q_e − Q_o)/4 and M = A^{-1/2}BA^{-1/2} (Cholesky whitening). Then A^{-1/2}(Q_e/2)A^{-1/2} = I + M and A^{-1/2}(Q_o/2)A^{-1/2} = I − M, so:
- 1 + μ_min(M) tests Q_e (the even sector);
- 1 − μ_max(M) tests Q_o (the odd sector).

This corrects the brief's pairing (it had them the other way round).

## Results
1. **Doubling identification (NUMERICAL, good to 3 digits).** λ_min(A(L)) = ½·min(q_e, q_o)(L/2) to 3 digits for ζ, DH and F_0.1 at L ≤ 1.1. For F_0.1 the counts also match exactly: n₋(A(L)) = n₋(Q(L/2)) = 1, 2, 2, 2, 3 at L 1.0–1.4. So A really is the half-window form. It drifts for ζ at L ≥ 1.2 (unconverged, below 1e−26) and for H₋ₜ at L ≥ 1.0. The H₋ₜ drift is probably zero truncation (MZ 50) or a span mismatch, and is flagged.
2. **Sign content (CHECKED, finite algebra).** ‖M(L)‖ ≤ 1 ⇔ Q_e(L), Q_o(L) ⪰ 0, given A(L) ≻ 0 ⇔ positivity at L/2. So the doubling step tells fakes apart exactly as window positivity at L does, and adds nothing. On this grid the first failures are: F_0.1 at L = 0.5 (q_e −3.9e−7, n₋ 1; PD at 0.45), H₋₄ at 0.5 (−1.1e−9), H₋₂ at 0.7 (−2.3e−15, LDLᵀ n₋ 1). DH stays PD throughout 0.3–1.4 (1 − ‖M‖ = 0.50 → 3.5e−14). ζ stays PD throughout, with 1 − ‖M‖ = 2.0e−2 → 4.0e−28 (L ≥ 1.0 is unconverged).
3. **Extremal M-eigenvector (CHECKED, expected).** For every PD form, the eigenvector at μ ≈ −1 is the Q_e near-null vector (cos ≥ 0.9999 for ζ; ≥ 0.98 for DH) and it is the S₊ near-null vector. The extremal side is always Q_e (μ → −1) for ζ and DH. Once a fake fails, the most negative generalised eigenvector moves away from the Q min-eigenvector (cos 0.01–0.8). That is a size or onset effect, not a sign separator.
4. **No doubling law.** 1 − ‖M(L)‖ is not fixed by half-window data:
   - (1 − ‖M‖)·λ_min(A)/(q_min(L)/2) = v_eᵀAv_e/λ_min(A) grows from 1.8 to 5e10 for ζ (L 0.3 → 1.4), 1.06 → 16 for DH, and 39–4e5 for F_0.1 and H₋₂ while they are still PD or near onset.
   - log(1 − ‖M(L)‖)/log(1 − ‖M(L/2)‖) is 4.8 → 2.3 for ζ and 2.6 → 10 for DH.
   - So there is no monotone or multiplicative law. Separation is by size only, and DH is never ruled out.

**Verdict: DEAD END.** H95a restates window positivity. The doubling or contraction view is the Kreĭn extension problem for positive-definite or screw functions (Kreĭn; Kreĭn–Langer for κ negative squares). Suzuki (arXiv 2206.03682, 2209.04658) has RH ⇔ the ζ screw function g ∈ 𝒢_∞, with positivity proved unconditionally on small intervals. Contraction-parametrised extensions are free in that theory, so a contraction bound needs an outside positivity source. CC 2106.01715 and 2310.18423 have no doubling statement.

Rules out: F_0.1 (L ≥ 0.5), H₋₄ (L ≥ 0.5), H₋₂ (L ≥ 0.7), all exactly where positivity does. DH is not ruled out.
