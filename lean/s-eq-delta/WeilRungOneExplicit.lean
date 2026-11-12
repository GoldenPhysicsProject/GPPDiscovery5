/-
Copyright (c) 2026 Daniel Toupin. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Daniel Toupin (with AI assistance: Heart of the Cosmos / Grok)
-/
import GppVerify.RiemannHypothesis.ArithmeticReflectionPositivity
import Mathlib.Analysis.Complex.ExponentialBounds

/-!
# Rung 1 of the Weil support ladder: the window, the explicit form, the certificate, the induction

Follow-up to `ArithmeticCFTRoutes.lean` (`arp_iff_all_rungs`: ARP ⇔ every rung). In-tree
convention: `ReflectionPositiveOnRung L` quantifies over smooth `f` supported in `[-L/2, L/2]`, so
`g = autocorr f` is supported in `[-L, L]` and the prime power `n` is visible iff `log n < L`.

## The window
* Rung 0: `L < log 2`, no prime power. Connes–Consani (cited; a named hypothesis in tree).
* **Rung 1: `log 2 ≤ L ≤ log 3`**, `f` supported in `[-log 3/2, log 3/2] ≈ [-0.5493, 0.5493]`. Only
  `n = 2` enters. (At `L = log 3` the point `g(log 3)` is the boundary value of a convolution square,
  so it is `0`.) Rung `k` ends at the `k+1`-st prime power: 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, …
* `weilForm_rung_one` (PROVED): on rung 1,
  `weilForm (autocorr f) = poleTerm g + archimedeanTerm g − (log 2/√2)·2·g(log 2)`,
  using `autocorr_neg` (PROVED: `g(−x) = g(x)`).

## The finite-dimensional problem
The rung-1 form is `Pole (rank 2) + A∞ (the log-Lévy Dirichlet form, ≥ C∞ = −5.372 and
growing like `log|t|`) − K₂`, where `K₂ = (log 2/√2)(τ_{log 2} + τ_{−log 2})` is a partial shift of norm ≤ `√2 log 2
≈ 0.980` on `L²[−h, h]`. Above the frequency `T₁` where `log(t/2π) − 1/t ≥ √2 log 2` the form is
positive, so the problem is a finite Galerkin matrix plus explicit tails (Zhu's one-stroke,
reproduced in GPPDiscovery5/s27).

## Status of rung 1 (labels)
* CERTIFIED (Arb ball arithmetic, single implementation, Lean-open; `s27_onestroke_arb.py` at
  `h = 11/20`, `T = 100`, `K = 60`): `Q(f) ≥ 3.479999999974342·10⁻⁸ ‖f‖²` (the recorded even `certified_lower`) for every `f` on `[−0.55, 0.55]`
  (even sector binding; odd `≥ 9.95999999999974·10⁻⁶`, recorded). Since `log 3 < 1.1` (`log_three_lt`, PROVED), this covers
  the whole closed rung 1. The earlier S27 certificate on `[−0.8, 0.8]` already covered it.
  `RungOneArbCertificate` is that statement as a named Prop. `rung_one_of_certificate` (PROVED) gives
  `ReflectionPositiveOnRung (log 3)`.
* NUMERICAL (float64, independent code, in-tree time-domain formulas): at the rung-1 edge
  `h = 0.5493`, `λ_min ≈ 5.58·10⁻⁸`, from a true-form Rayleigh quotient of `5.5807·10⁻⁸` against
  `λ_min(R_T) = 5.5715·10⁻⁸` at `T = 3200`, which increases with `T`. The extremal test is even,
  nodeless and bell-shaped (Perron–Frobenius ground state). Budget: `pole/‖f‖² = 1.3044`,
  `A∞/‖f‖² = −1.2865`, `kick/‖f‖² = 0.0178` (`g(log 2)/g(0) = 0.0182`, against the Cauchy–Schwarz
  bound 1). The kick is about 3·10⁵ times the margin, and the extremal test keeps `g(log 2)` tiny.
  So rung 1 is **not** a perturbation of rung 0: the sufficient condition `pole + A∞ ≥ √2 log 2 ·
  g(0)` is false on the extremal.
* Lean proof of rung 1: OPEN. It would need the digamma envelope, quadrature error bounds and a
  verified PD check. Rung 0 itself is not in Lean either.

## Margin across rungs (lower bounds CERTIFIED by Arb; the rest NUMERICAL)
| rung edge | window `L` | certified `λ_min ≥` | NUMERICAL `λ_min ≈` |
|---|---|---|---|
| 0 | `log 2` | 1.2999999999999997e−3 | 1.331e−3 |
| 1 | `log 3` | 4.8099999999746e−8 (h = 0.5493) | 5.58e−8 |
| 2 | `log 4` | 5.79545426360329e−13 | 8.5e−13 (float floor) |
| 3 (S27) | 1.6 < `log 5` | 1.02e−17 | — |
| (S27/28) | 2.0 | 4.12e−30 | — |
The fit is `ln λ_min ≈ −c·e^L + d` (NUMERICAL, from the three NUMERICAL margins above): successive
slopes `−10.08` (rung 0→1) and `−11.09` (rung 1→2), least-squares `c ≈ 10.6`, `d ≈ +14.7`. So `c ≈ 11`
is the slope, not a pure `exp(−11·e^L)` (that would give `≈ 5·10⁻¹⁵` at rung 1, not `5.58·10⁻⁸`).
The margin is doubly exponential in the window. Certified bounds are stated at the recorded
`certified_lower` digits (Turbo 6815936710, correction 1).

## Induction k → k+1 and where it breaks
* The form does not change between rungs; only the test class grows. So rung `k+1` ⇒ rung `k`
  (`rung_mono`, PROVED), and nothing propagates upward.
* Splitting tests as (old window) ⊕ (new directions), the step is a Schur complement. In the 2×2
  model (`two_by_two_nonneg_iff`, PROVED), with `a = λ_k` the rung-`k` margin, `d` the diagonal of a
  new direction and `b` its coupling to the soft mode, positivity holds iff `b² ≤ λ_k·d`. With
  `ln λ_k ≈ −c e^{L_k} + d` the coupling of every new direction to the rung-`k` ground state must be
  doubly-exponentially small. Nothing in rung `k` supplies that; it is a fresh arithmetic
  cancellation at every rung. **This is where the induction breaks**: rung `k` gives no quantitative
  input for rung `k+1`. `RungStepCoupling` names the needed bound.

Labels: theorems PROVED (standard axioms). `RungOneArbCertificate` is CERTIFIED externally (Arb)
and named here, not asserted. `RungStepCoupling` and Lean rung 1 are OPEN. ARP OPEN. RH OPEN.

Daniel Toupin / Golden Physics Project … This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
Audit: Turbo.
-/

open MeasureTheory
open scoped ContDiff
open GppArithmeticReflection

noncomputable section

namespace GppWeilRungOneExplicit

/-- **The autocorrelation is even:** `g(−x) = g(x)`. PROVED. -/
theorem autocorr_neg (f : ℝ → ℝ) (x : ℝ) : autocorr f (-x) = autocorr f x := by
  unfold autocorr
  have h := integral_sub_right_eq_self (μ := (volume : Measure ℝ))
    (fun y => f y * f (y + x)) x
  simp only [sub_neg_eq_add, sub_add_cancel] at h ⊢
  rw [← h]
  congr 1; funext y; ring

/-- **The explicit rung-1 form.** For `f` supported in `[-L/2, L/2]` with `L < log 3`, only `n = 2`
is visible: `weilForm g = pole g + A∞ g − (log 2/√2)·2·g(log 2)`. PROVED. -/
theorem weilForm_rung_one {f : ℝ → ℝ} {L : ℝ} (hsupp : ∀ y : ℝ, L / 2 < |y| → f y = 0)
    (hL : L < Real.log 3) :
    weilForm (autocorr f) = poleTerm (autocorr f) + GppArchimedeanFloor.archimedeanTerm (autocorr f) -
      Real.log 2 / Real.sqrt 2 * (2 * autocorr f (Real.log 2)) := by
  have hN : L < Real.log ((2 : ℕ) + 1) := by norm_num; exact hL
  rw [weilForm_eq_truncation hsupp hN]
  simp only [Finset.sum_range_succ, Finset.range_zero, Finset.sum_empty,
    ArithmeticFunction.vonMangoldt_apply_prime Nat.prime_two]
  have h0 : ArithmeticFunction.vonMangoldt 0 = 0 := by simp
  have h1 : ArithmeticFunction.vonMangoldt 1 = 0 := by simp
  rw [h0, h1]
  push_cast
  rw [autocorr_neg f (Real.log 2)]
  ring

/-- **Monotonicity of rungs:** the form is the same at every rung; only the test class grows.
PROVED. -/
theorem rung_mono {L L' : ℝ} (h : ReflectionPositiveOnRung L) (hle : L' ≤ L) :
    ReflectionPositiveOnRung L' := fun f hf hc hs =>
  h f hf hc fun y hy => hs y (by linarith)

/-- `log 3 < 11/10` (so the certificate window `[−0.55, 0.55]` covers rung 1). PROVED. -/
theorem log_three_lt : Real.log 3 < 11 / 10 := by
  rw [Real.log_lt_iff_lt_exp (by norm_num)]
  have h1 := Real.exp_one_gt_d9
  have h2 := Real.quadratic_le_exp_of_nonneg (x := 1 / 10) (by norm_num)
  have : Real.exp (11 / 10) = Real.exp 1 * Real.exp (1 / 10) := by
    rw [← Real.exp_add]; norm_num
  rw [this]
  nlinarith [Real.exp_pos (1 / 10)]

/-- **The rung-1 Arb certificate, as a named Prop** (CERTIFIED externally by ball arithmetic,
`s27_onestroke_arb.py --L 11/20 --T 100 --K 60`; single implementation; not a Lean proof):
`weilForm (autocorr f) ≥ 3.47·10⁻⁸ ∫ f²` for every smooth `f` supported in `[−11/20, 11/20]`.
The recorded `certified_lower` is `3.479999999974342·10⁻⁸`, just below `3.48·10⁻⁸`, so the Prop uses
`347/10¹⁰`, which the artifact certifies exactly (Turbo 6815936710, correction 1). -/
def RungOneArbCertificate : Prop :=
  ∀ f : ℝ → ℝ, ContDiff ℝ ∞ f → HasCompactSupport f → (∀ y : ℝ, 11 / 20 < |y| → f y = 0) →
    (347 / 10 ^ 10 : ℝ) * ∫ y, f y ^ 2 ≤ weilForm (autocorr f)

/-- **Certificate ⇒ positivity on the whole closed rung 1** (`L = log 3`). PROVED (conditional on
the external certificate). -/
theorem rung_one_of_certificate (h : RungOneArbCertificate) :
    ReflectionPositiveOnRung (Real.log 3) := by
  have hwide : ReflectionPositiveOnRung (11 / 10) := fun f hf hc hs => by
    have := h f hf hc fun y hy => hs y (by linarith)
    have hint : 0 ≤ ∫ y, f y ^ 2 := integral_nonneg fun y => sq_nonneg _
    nlinarith
  exact rung_mono hwide log_three_lt.le

/-- **Schur step in the 2×2 model:** for `a > 0`, `a x² + 2bxy + d y² ≥ 0` for all `x, y` iff
`b² ≤ a d`. PROVED. -/
theorem two_by_two_nonneg_iff {a b d : ℝ} (ha : 0 < a) :
    (∀ x y : ℝ, 0 ≤ a * x ^ 2 + 2 * b * x * y + d * y ^ 2) ↔ b ^ 2 ≤ a * d := by
  constructor
  · intro h
    have := h (-b) a
    have h' : 0 ≤ a * (a * d - b ^ 2) := by nlinarith
    nlinarith [(mul_nonneg_iff_of_pos_left ha).1 h']
  · intro hb x y
    have key : a * (a * x ^ 2 + 2 * b * x * y + d * y ^ 2) =
        (a * x + b * y) ^ 2 + (a * d - b ^ 2) * y ^ 2 := by ring
    have : 0 ≤ a * (a * x ^ 2 + 2 * b * x * y + d * y ^ 2) := by
      rw [key]; nlinarith [sq_nonneg (a * x + b * y), sq_nonneg y]
    exact (mul_nonneg_iff_of_pos_left ha).1 this

/-- **The induction step needs a coupling bound** (named, OPEN). For the rung-`k` margin `lamK`,
every new direction's diagonal `d` and coupling `b` to the rung-`k` soft mode satisfy `b² ≤ lamK·d`.
In the 2×2 model this is exactly the step (`two_by_two_nonneg_iff`). -/
def RungStepCoupling (lamK : ℝ) (couplings : Set (ℝ × ℝ)) : Prop :=
  ∀ bd ∈ couplings, bd.1 ^ 2 ≤ lamK * bd.2

/-- In the model, the coupling bound is equivalent to positivity of every 2×2 step. PROVED. -/
theorem rungStepCoupling_iff {lamK : ℝ} (hl : 0 < lamK) (couplings : Set (ℝ × ℝ)) :
    RungStepCoupling lamK couplings ↔
      ∀ bd ∈ couplings, ∀ x y : ℝ, 0 ≤ lamK * x ^ 2 + 2 * bd.1 * x * y + bd.2 * y ^ 2 :=
  ⟨fun h bd hbd => (two_by_two_nonneg_iff hl).2 (h bd hbd),
    fun h bd hbd => (two_by_two_nonneg_iff hl).1 (h bd hbd)⟩

end GppWeilRungOneExplicit

#print axioms GppWeilRungOneExplicit.autocorr_neg
#print axioms GppWeilRungOneExplicit.weilForm_rung_one
#print axioms GppWeilRungOneExplicit.rung_mono
#print axioms GppWeilRungOneExplicit.log_three_lt
#print axioms GppWeilRungOneExplicit.rung_one_of_certificate
#print axioms GppWeilRungOneExplicit.two_by_two_nonneg_iff
#print axioms GppWeilRungOneExplicit.rungStepCoupling_iff
