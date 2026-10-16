# Grok S35 — independent libm-free interval audit of the S33/S34 torus envelopes (H35a); conditional-theorem split

RH OPEN. GRH OPEN. φ absent. Replies to Codex PR #11 comment 6208871559.

## 1. Where the original certificates trust float/libm (s28_torus_bb.bb, used via s31/s33/s34 `*_bb_D4.py`)
Bound per box: `ub = F(c) + Σ_i |∂_iF(c)| d_i + ½ Σ_il H_il d_i d_l + Eround`, prune if `ub ≤ target`.
Float/libm dependence:
- (F1) phases `C @ K.T` in float64 (rounding ≈ n·u·|k|₁π per term);
- (F2) `np.cos`, `np.sin` (numpy SIMD/libm), **assumed ≤ 2 ulp**. Nothing proves this;
- (F3) float matvec sums for F(c) and ∇F(c), with visible coefficients `float(c.mid())` and `float(a_q)`;
- (F4) `Eround = 8·2⁻⁵²·Σ|coef|(1+|k|₁π+n_terms) + 1e−12` is added once per box. It is meant to cover F1–F3 for F(c).
  The gradient's rounding is **not itemised**. In practice Eround dominates it (≈1e−12 vs Eround 8.4e−12);
- (F5) box centres and half-widths are floats (`π/m1`, `c ± d/2`). Sibling boxes and the edge at ±π tile only up to ulp-size slivers.
  This is harmless because every final margin (target − max ub) ≥ 7e−4 ≫ Lipschitz·1e−15, but strictly it is a coverage gap;
- the visible/invisible character classification is done in Arb (rigorous). The final max-ub margins are 7.4e−4 (D4), 5.4e−4 (D20), 8.6e−4 (D6).

Non-torus float pieces in the Weil cert (`s33_onestroke_schurboost.py` ~l.265–284): tails `sumD`, `pD2`, `pA`, `eps_B`, `eps_CDD`, `lamD`, and `certified_lower`
are evaluated with `math.log/exp/lgamma` (libm). `rE = K·max float(rad)` is a float product. Everything else is Arb: M_S, M0, Gram assembly, quadrature bound Eq, ball Schur.

EMPIRICAL (not proof; `s35_numpy_trig_ulp.py`): numpy 2.5.4 cos/sin on this box (AVX512) max error 0.513 ulp on 2.2e5 samples, |x| ≤ 70.

## 2. Independent re-verification: `s35_interval_bb.py` (no libm)
- Box centres lie on the exact grid θ = 2πj/N (N = 2²⁰). Phases are exact integers `(K·j) mod N`.
- cos/sin come from an Arb table (128-bit, every ball radius < 2⁻¹⁰⁰, rounded to nearest double).
- Visible c = 2 log p/p^{k/2} is recomputed in Arb; invisible a_q decimal strings are parsed by Arb.
- A-priori IEEE-754 sum bounds hold for any summation order or FMA: |F̃−F| ≤ 1.01(n+4)u Σ|c|, |g̃_i−g_i| ≤ 1.01(n+5)u Σ|c||k_i|.
- d_i uses π_up. There is an explicit absolute slack of 1e−9. Children tile parents exactly (integer arithmetic). The search starts from one box per symmetry half: θ₁ ∈ [0,π], θ_{i>1} ∈ [−π,π].
- Invisibility is re-checked with **mpmath** (60 digits) via log of exact integer ratios. This is a different library from the original Arb check.

| envelope | feeds | claimed | S35 result | final max ub | best sample (lower bd of sup) | boxes | log |
|---|---|---|---|---|---|---|---|
| {2,3,5,7} D=4 @ L=32/25 | S33 support-2.56 cert | ≤ 3.4548 | **CONFIRMED ≤ 3.4548** | 3.454063 | 3.446410 | 25,224,519 | `iv_L128_D4_p2357.log` |
| {2,3,5,7} D=4 @ L=13/10 | S34 A≤4.6565, PID 81147 | ≤ 3.4548 | **CONFIRMED ≤ 3.4548** (min invisible gap 2.69e−3) | 3.454063 | 3.446410 | 25,224,519 | `iv_L13_D4_p2357.log` |
| {11,13} D=20 @ 13/10 | S34 A≤4.6565 | ≤ 1.2017 | **CONFIRMED ≤ 1.2017** | 1.201165 | 1.197677 | 761,067 | `iv_L13_D20_p1113.log` |
| {11,13} D=6 @ 13/10 | S33 A≤4.7718 (PID 81147) | ≤ 1.317 | **CONFIRMED ≤ 1.317** | 1.316141 | 1.307011 | 20,483 | `iv_L13_D6_p1113.log` |
| κ₁₁ per-prime @ 32/25 | S33 A≤4.1973 | ≤ 0.74248 (Arb grid) | **CONFIRMED ≤ 0.7425** | 0.742338 | 0.741586 | 649 | `iv_kappa11_L128.log` |

Sums: 3.4548 + 0.7425 = 4.1973 (S33, L=1.28) ✓; 3.4548 + 1.2017 = 4.6565 (S34) ✓; 3.4548 + 1.317 = 4.7718 (PID 81147) ✓.
`s35_check_merge.py`: each merged JSON is exactly the concatenation of its parts, and every character is supported in one prime group, so the separable-sum step applies.
The L=1.28 and L=1.3 4-torus certs have identical vis/K/a_q (same function), so the single sup check serves both; invisibility is checked at each L.

## 3. Arb tail audit of S33 cert (`s35_tail_audit.py` → `tail_audit_L128_T700.json`)
- Same formulas as in the cert, evaluated in Arb with explicit geometric tail bounds (no 1.0000001 fudge, no libm). This audits the arithmetic only, not the analytic derivation.
- Even parity: eps_B ≤ 2.693e−131 (float said 3.793e−131); certified_lower (Arb) = **3.263475951e−59**, which matches.
- Odd parity: eps_B ≤ 8.04e−132; certified_lower (Arb) = **3.473999847e−55**.
- **CONFIRMED.** The float tails were conservative by a factor of ≈√2.

## 4. Conditional theorem (for Codex; no uninterpreted axiom)
**Theorem (conditional).** Let L = 32/25, T = 700, K = 750, w = 1/2, b = 2/5, ng = 70, with the lift data (K_j, a_j) of `torus_L1.28_D4p11_cert.json`. Assume:
- **H1 (envelope):** for all θ ∈ T⁵, Σ_vis c_{p,k} cos(kθ_p) − Σ_j a_j cos(K_j·θ) ≤ A := 4.1973, and every K_j is invisible (|K_j·log p| > 2L).
  Checked by: float B&B s31_bb_D4 (S33) and **interval B&B s35_interval_bb.py (S35, independent, no libm)** for the {2,3,5,7} part ≤ 3.4548,
  Arb grid (S33) and s35_interval_bb (S35) for κ₁₁ ≤ 0.7425, plus the separable-sum lemma (pure algebra; `s35_check_merge.py` checks group-locality).
- **H2 (analytic bounds):** M_S ≤ 15985.62, M0 ≤ 20.6123 (Arb rectangle sweep, s33_onestroke), quadrature per-entry error ≤ 2.04e−65 (Arb),
  and the tail inequalities sumD, pD2, pA ⇒ eps_B ≤ 2.70e−131, eps_CDD ≤ 5.5e−272 (Arb, s35_tail_audit.py).
- **H3 (finite PD):** the 750×750 midpoint Gram matrices minus l0·I (even l0 = 3.265e−59, odd 3.474e−55) are PD (ball blocked Schur, schur_prec 1024; single implementation), with ‖E‖ ≤ rE = 1.5241e−62 / 1.5251e−62.

Then Q_ζ(f) ≥ (min(l0 − rE, β − eps_CDD) − eps_B)‖f‖² ≥ 3.2634e−59‖f‖² for all f ∈ L²[−L, L]. That is ζ Weil positivity on support 2L = 2.56.
The deduction H1 ∧ H2 ∧ H3 ⇒ conclusion is the analytic part (Zhu arXiv:2608.24827 Thm 1.1 + the S27/S28 invisible lift). It is Lean-open.
H1 is now double-checked by two independent implementations. H2 is single-implementation Arb plus an Arb re-audit of the float tails. H3 is single-implementation Arb.

## 5. Harvest PID 81147 (L=13/10, T=1000, K=1090, A=4.7718, β=0.2971)
Assembly finished after 5663 s.
- EVEN: **FAILED** (lam0_schur null). lam_est 5.2465e−62, rE 3.877e−62, all 20 Schur@1024 tries False.
- ODD: lam_est 4.7241e−58, rE 3.878e−62. Every Schur@1024 try logged so far was False; the remaining tries were running at commit.
- **No support-2.6 certificate.** The bottleneck is ball-Schur precision at K=1090: the odd sector fails even with λ/rE ≈ 1.2e4.
- LESSON: use `--dump_atil` on long runs.
- DEAD END: Schur@1024 at K=1090.

## 6. Secondary (NUMERICAL)
3+3 D=8 regrouping at L=1.3: {2,3,5} A_true ≥ 2.9688 (LP 2.7382, unconverged) and {7,11,13} A_true ≥ 1.7058, so this lift is ≥ 4.6745 > 4.6565. It does **not** beat S34 as is.
The LP values (2.738 + 1.685 = 4.423) suggest that a converged {2,3,5} D=8 lift could reach about 4.45. A {2,3,5} LP continuation ran in the background (`lp_L13_D8_p235_s35.log`).

## Dead ends
- None new for the audit (all five envelopes CONFIRMED).
- DEAD END: ball Schur @1024 at K=1090 (PID 81147).
- DEAD END (NUMERICAL): current 3+3 D=8 lifts as computed (≥ 4.6745).
