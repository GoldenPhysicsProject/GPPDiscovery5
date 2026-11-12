/-
Copyright (c) 2026 Daniel Toupin. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Daniel Toupin (with AI assistance: Heart of the Cosmos / Grok)
-/
import GppVerify.RiemannHypothesis.ArithmeticReflectionPositivity
import Mathlib.Analysis.Complex.ExponentialBounds

/-!
# Rung 2 of the Weil support ladder and `RungStepCoupling` written out (n = 3, then n = 4)

Follow-up to `WeilRungOneExplicit.lean` (GPPDiscovery5 `cc284e4`, corrected in `f7bed28`).
In-tree convention: `ReflectionPositiveOnRung L` means `f` is supported in `[-L/2, L/2]`.

## Windows and explicit forms (PROVED)
* Rung 2 is `log 3 ≤ L ≤ log 4`, i.e. `f` on `[-log 2, log 2]`; `n = 2, 3` are visible.
  `weilForm_rung_two`: `W = pole + A∞ − (log 2/√2)·2g(log 2) − (log 3/√3)·2g(log 3)`.
* Rung 3 begins at `L = log 4`, where `n = 4` enters (`Λ(4) = log 2`, `√4 = 2`).
  `weilForm_rung_three`: for `L < log 5`, add `−(log 2/2)·2g(log 4)`.
* `log_four_lt`: `log 4 < 7/5`, so a certificate on `[-7/10, 7/10]` covers the closed rung 2.

## The step as a block form (PROVED, abstract)
Split the new window as `V_k ⊕ V_k^⊥`, where `V_k^⊥` is supported on the annulus `h_k < |x| ≤ h_{k+1}`.
Let `B` be the polarization of `f ↦ weilForm (autocorr f)`. (The bilinear structure is STANDARD
but is not formalized here, so the block lemmas below are stated for an abstract symmetric
bilinear form on a real inner product space.)
* `sq_le_mul_of_forall_nonneg`: if `a x² + 2bxy + d y² ≥ 0` for all `x, y`, then `b² ≤ a d`.
  No sign hypothesis is needed.
* `coupling_le_of_margin`: if `B v v ≥ μ‖v‖²` on the whole new window, then for orthonormal
  `φ ∈ V_k`, `ψ ∈ V_k^⊥`: `B(φ,ψ)² ≤ (B φφ − μ)(B ψψ − μ)`.
* `rungStepCoupling_of_margin`: if moreover `φ` is the rung-`k` ground state, `B φφ = λ_k`, then
  `b² ≤ λ_k·d` (`RungStepCoupling`) holds for every unit `ψ ⊥ φ`.
  This runs **downward**: rung `k+1` ⇒ the coupling bound. The coupling bound is the 2×2
  necessary condition for the upward step, not a sufficient one (the full step is the Schur
  complement `D − Bᵀ A⁻¹ B ≥ 0`).

## The numbers, step rung 1 → 2 (`n = 3` enters). No zeros as input.
CERTIFIED (Arb ball arithmetic, single implementation, Lean-open; bounds quoted at their recorded
digits, never rounded up):
* `λ₁ ≥ 3.479999999974342e-8`: the rung-1 certificate on `[-11/20, 11/20] ⊇ [-log 3/2, log 3/2]`.
* **Closed rung 2:** for `f` on `[-7/10, 7/10] ⊇ [-log 2, log 2]` (comb `{2,3,4}`, `T = 150`,
  `K = 90`, prec 256), `Q(f) ≥ 2.452756856955386e-13 ‖f‖²`. Even sector binding; odd
  `≥ 1.4499920447324892e-10`. The earlier `h = 6931/10000` certificate (`5.79545426360329e-13`)
  stops just short of `log 2 = 0.693147…`, so it does not cover the closed rung. A
  `T = 200`, `K = 110` run failed in the even sector (`certified_lower < 0`, tail coupling `eps_B`
  too large) and is recorded but not used. In Lean the certificate is the named Prop
  `RungTwoArbCertificate` (constant `245/10^15`, below the recorded bound), and
  `rung_two_of_certificate` (PROVED) turns it into `ReflectionPositiveOnRung (log 4)`.
* Through `coupling_le_of_margin`, the rung-2 certificate implies the coupling bound for every
  orthonormal pair with margin `μ = 2.45·10⁻¹³`. This is CERTIFIED only as a **consequence** of the
  direct rung-2 certificate. There is no independent (upward) proof of it.

NUMERICAL (float64, `num/rung2/rung_coupling.py`; Galerkin of `R_T` at `T = 800, 3200`, with
Legendre bases of 44 modes on `V₁` and 16 on each annulus piece; true-form values in the time domain
from the in-tree definitions):
* `λ₁ ≈ 5.58e-8` (`Q(φ₁) = 5.5821e-8`). `φ₁` is even, nodeless and bell-shaped.
* Worst new direction for the 2×2 test, `ψ_w ∝ D⁻¹ Bᵀ φ₁` on the annulus `log 3/2 < |x| ≤ log 2`:
  `b = Q(φ₁, ψ_w) ≈ 2.99e-4` and `d = Q(ψ_w, ψ_w) ≈ 2.177`, so
  **`b²/(λ₁ d) ≈ 0.735 < 1`** (`R_3200`: 0.7350; true form: 0.7354). So `RungStepCoupling` holds with
  about 26% room for every new direction.
* Split of `b` (true form): pole `−0.21730`, A∞ `+0.09714`, n=2 kick `+0.11086`, n=3 kick
  `+0.00960`. That is a cancellation by a factor of about 700. **The n=3 kick is about 32·b.**
  Without the `n = 3` term, `b ≈ −0.0093` and `b²/(λ₁ d) ≈ 700`, so the 2×2 test fails by
  orders of magnitude. The new prime is exactly what makes the step possible.
* Schur complement on the annulus: `λ_min(D − Bᵀ A⁻¹ B) ≈ 1.69e-7 > 0`. Block minimum:
  `λ₂ ≈ 7.5e-13` (Galerkin, not converged in the last digit). The certified lower bound is above.

## The numbers, step rung 2 → 3 (`n = 4` enters). NUMERICAL at about 1 digit
* `λ₂ ≈ 7.6e-13 … 8.9e-13` (it depends on the basis; this is near the float64 floor).
* Worst direction on the annulus `log 2 < |x| ≤ log 5/2`: `b ≈ 1.1e-6 … 1.2e-6`, `d ≈ 2.9`, so
  `b²/(λ₂ d) ≈ 0.54 … 0.58 < 1`.
* Split of `b`: pole `+0.0586`, A∞ `−0.0234`, n=2 `−0.0215`, n=3 `−0.0136`, n=4 `−7.7e-5`.
  **The n=4 kick is about 68·b.** Without it, `b²/(λ₂ d) ≈ 2–3·10³`.
* `λ₃` and the Schur complement lie below the float64 floor: the float run returns
  `λ_min ≈ −1.9e-14` and Schur `≈ −3e-3`, which are artefacts of `A⁻¹` with condition about `10¹²`.
  NOT resolved here. Rung 3 rests on S27 Arb certificates. The `h = 0.8` certificate
  (`1.02e-17`) covers `L ≤ 1.6 < log 5`, so it does not reach the end of rung 3. The `h = 1`
  certificate (`4.13e-30`) covers all of rung 3. `num/rung1/cert_rung3_edge.json` has
  `certified_lower = −1894.5`: it certifies nothing.

## What this says about the induction
At both steps the coupling `b` is the residue of a cancellation among the pole, A∞ and the old
primes, and the new prime's kick is 30–70 times larger than `b`. The 2×2 test passes with
`O(1)` room (0.74, then about 0.56). But `b` itself must be `O(√λ_k)`, i.e. doubly exponentially
small. Nothing at rung `k` predicts the size of the new prime's contribution, so the step is a
fresh arithmetic identity at each rung (`RungStepCoupling`, OPEN in general; true at the n=3 step
as a consequence of the rung-2 certificate).

Labels: theorems PROVED (standard axioms). `RungTwoArbCertificate` is CERTIFIED externally (Arb)
and named, not asserted. `RungStepCoupling` is OPEN in general. ARP OPEN. RH OPEN.

Daniel Toupin / Golden Physics Project … This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
Audit: Turbo.
-/

open MeasureTheory
open scoped ContDiff RealInnerProductSpace
open GppArithmeticReflection

noncomputable section

namespace GppWeilRungTwo

/-- The autocorrelation is even (local copy of `GppWeilRungOneExplicit.autocorr_neg`). PROVED. -/
theorem autocorr_neg (f : ℝ → ℝ) (x : ℝ) : autocorr f (-x) = autocorr f x := by
  unfold autocorr
  have h := integral_sub_right_eq_self (μ := (volume : Measure ℝ))
    (fun y => f y * f (y + x)) x
  simp only [sub_neg_eq_add, sub_add_cancel] at h ⊢
  rw [← h]
  congr 1; funext y; ring

/-- **Explicit rung-2 form:** for `L < log 4`, only `n = 2, 3` are visible. PROVED. -/
theorem weilForm_rung_two {f : ℝ → ℝ} {L : ℝ} (hsupp : ∀ y : ℝ, L / 2 < |y| → f y = 0)
    (hL : L < Real.log 4) :
    weilForm (autocorr f) = poleTerm (autocorr f) + GppArchimedeanFloor.archimedeanTerm (autocorr f) -
      Real.log 2 / Real.sqrt 2 * (2 * autocorr f (Real.log 2)) -
      Real.log 3 / Real.sqrt 3 * (2 * autocorr f (Real.log 3)) := by
  have hN : L < Real.log ((3 : ℕ) + 1) := by norm_num; exact hL
  rw [weilForm_eq_truncation hsupp hN]
  simp only [Finset.sum_range_succ, Finset.range_zero, Finset.sum_empty,
    ArithmeticFunction.vonMangoldt_apply_prime Nat.prime_two,
    ArithmeticFunction.vonMangoldt_apply_prime Nat.prime_three]
  have h0 : ArithmeticFunction.vonMangoldt 0 = 0 := by simp
  have h1 : ArithmeticFunction.vonMangoldt 1 = 0 := by simp
  rw [h0, h1]
  push_cast
  rw [autocorr_neg f (Real.log 2), autocorr_neg f (Real.log 3)]
  ring

/-- **Explicit form at the start of rung 3:** for `L < log 5`, `n = 4` (`Λ(4) = log 2`) joins.
PROVED. -/
theorem weilForm_rung_three {f : ℝ → ℝ} {L : ℝ} (hsupp : ∀ y : ℝ, L / 2 < |y| → f y = 0)
    (hL : L < Real.log 5) :
    weilForm (autocorr f) = poleTerm (autocorr f) + GppArchimedeanFloor.archimedeanTerm (autocorr f) -
      Real.log 2 / Real.sqrt 2 * (2 * autocorr f (Real.log 2)) -
      Real.log 3 / Real.sqrt 3 * (2 * autocorr f (Real.log 3)) -
      Real.log 2 / 2 * (2 * autocorr f (Real.log 4)) := by
  have hN : L < Real.log ((4 : ℕ) + 1) := by norm_num; exact hL
  rw [weilForm_eq_truncation hsupp hN]
  have h4 : ArithmeticFunction.vonMangoldt 4 = Real.log 2 := by
    have : (4 : ℕ) = 2 ^ 2 := by norm_num
    rw [this, ArithmeticFunction.vonMangoldt_apply_pow (by norm_num),
      ArithmeticFunction.vonMangoldt_apply_prime Nat.prime_two]
    simp
  simp only [Finset.sum_range_succ, Finset.range_zero, Finset.sum_empty,
    ArithmeticFunction.vonMangoldt_apply_prime Nat.prime_two,
    ArithmeticFunction.vonMangoldt_apply_prime Nat.prime_three, h4]
  have h0 : ArithmeticFunction.vonMangoldt 0 = 0 := by simp
  have h1 : ArithmeticFunction.vonMangoldt 1 = 0 := by simp
  have hs4 : Real.sqrt (4 : ℕ) = 2 := by
    rw [show ((4 : ℕ) : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]
  rw [h0, h1, hs4]
  push_cast
  rw [autocorr_neg f (Real.log 2), autocorr_neg f (Real.log 3), autocorr_neg f (Real.log 4)]
  ring

/-- `log 4 < 7/5`. PROVED. -/
theorem log_four_lt : Real.log 4 < 7 / 5 := by
  have h : Real.log 4 = 2 * Real.log 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; norm_num
  rw [h]; linarith [Real.log_two_lt_d9]

/-- **The rung-2 Arb certificate, as a named Prop** (CERTIFIED externally:
`s27_onestroke_arb.py --L 7/10 --T 150 --K 90 --ng 24 --prec 256`, recorded even
`certified_lower = 2.452756856955386e-13`, odd `1.4499920447324892e-10`; single implementation;
not a Lean proof). -/
def RungTwoArbCertificate : Prop :=
  ∀ f : ℝ → ℝ, ContDiff ℝ ∞ f → HasCompactSupport f → (∀ y : ℝ, 7 / 10 < |y| → f y = 0) →
    (245 / 10 ^ 15 : ℝ) * ∫ y, f y ^ 2 ≤ weilForm (autocorr f)

/-- **Certificate ⇒ positivity on the whole closed rung 2** (`L = log 4`). PROVED (conditional). -/
theorem rung_two_of_certificate (h : RungTwoArbCertificate) :
    ReflectionPositiveOnRung (Real.log 4) := by
  have hwide : ReflectionPositiveOnRung (7 / 5) := fun f hf hc hs => by
    have := h f hf hc fun y hy => hs y (by linarith)
    have hint : 0 ≤ ∫ y, f y ^ 2 := integral_nonneg fun y => sq_nonneg _
    nlinarith
  exact fun f hf hc hs => hwide f hf hc fun y hy => hs y (by linarith [log_four_lt])

/-- **2×2 Cauchy–Schwarz without sign hypotheses:** a nonnegative binary quadratic form has
`b² ≤ a d`. PROVED. -/
theorem sq_le_mul_of_forall_nonneg {a b d : ℝ}
    (h : ∀ x y : ℝ, 0 ≤ a * x ^ 2 + 2 * b * x * y + d * y ^ 2) : b ^ 2 ≤ a * d := by
  have ha : 0 ≤ a := by simpa using h 1 0
  rcases ha.lt_or_eq with ha | ha
  · have := h (-b) a
    have h' : 0 ≤ a * (a * d - b ^ 2) := by nlinarith
    nlinarith [(mul_nonneg_iff_of_pos_left ha).1 h']
  · subst ha
    by_contra hb
    have hb0 : b ≠ 0 := by rintro rfl; simp at hb
    have := h (-(d + 1) / (2 * b)) 1
    have e : 2 * b * (-(d + 1) / (2 * b)) * 1 = -(d + 1) := by field_simp
    linarith [this, e]

/-- **2×2 margin lemma:** if the binary form is `≥ μ(x² + y²)`, then `b² ≤ (a − μ)(d − μ)`.
PROVED. -/
theorem sq_le_of_margin {a b d μ : ℝ}
    (h : ∀ x y : ℝ, μ * (x ^ 2 + y ^ 2) ≤ a * x ^ 2 + 2 * b * x * y + d * y ^ 2) :
    b ^ 2 ≤ (a - μ) * (d - μ) :=
  sq_le_mul_of_forall_nonneg fun x y => by nlinarith [h x y]

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- **Coupling bound from a window margin.** If `B v v ≥ μ‖v‖²` on the whole new window, then for
orthonormal `φ` (old window) and `ψ` (annulus): `B(φ,ψ)² ≤ (B φφ − μ)(B ψψ − μ)`. PROVED. -/
theorem coupling_le_of_margin (B : LinearMap.BilinForm ℝ E) (hB : ∀ u v, B u v = B v u) {μ : ℝ}
    (hμ : ∀ v, μ * ‖v‖ ^ 2 ≤ B v v) {φ ψ : E} (hφ : ‖φ‖ = 1) (hψ : ‖ψ‖ = 1)
    (horth : ⟪φ, ψ⟫ = 0) :
    (B φ ψ) ^ 2 ≤ (B φ φ - μ) * (B ψ ψ - μ) := by
  apply sq_le_of_margin
  intro x y
  have hn : ‖x • φ + y • ψ‖ ^ 2 = x ^ 2 + y ^ 2 := by
    rw [← real_inner_self_eq_norm_sq]
    simp only [inner_add_left, inner_add_right, real_inner_smul_left, real_inner_smul_right,
      real_inner_self_eq_norm_sq, horth, real_inner_comm φ ψ]
    rw [norm_smul, norm_smul, Real.norm_eq_abs, Real.norm_eq_abs, mul_pow, mul_pow, sq_abs, sq_abs,
      hφ, hψ]
    ring
  have hq : B (x • φ + y • ψ) (x • φ + y • ψ) =
      B φ φ * x ^ 2 + 2 * B φ ψ * x * y + B ψ ψ * y ^ 2 := by
    simp only [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul, hB ψ φ]
    ring
  have := hμ (x • φ + y • ψ)
  rw [hn, hq] at this
  linarith

/-- **`RungStepCoupling` (named, as in `WeilRungOneExplicit`):** every listed (coupling,
diagonal) pair satisfies `b² ≤ λ_k·d`. OPEN in general. -/
def RungStepCoupling (lamK : ℝ) (couplings : Set (ℝ × ℝ)) : Prop :=
  ∀ bd ∈ couplings, bd.1 ^ 2 ≤ lamK * bd.2

/-- **The coupling bound at a step follows (downward) from a margin on the bigger window.** If
`B ≥ μ ≥ 0` on the new window and `φ` is a unit rung-`k` mode with `B φφ = λ_k`, then
`RungStepCoupling λ_k` holds for the couplings of every unit `ψ ⊥ φ`. PROVED. -/
theorem rungStepCoupling_of_margin (B : LinearMap.BilinForm ℝ E) (hB : ∀ u v, B u v = B v u)
    {μ : ℝ} (hμ0 : 0 ≤ μ) (hμ : ∀ v, μ * ‖v‖ ^ 2 ≤ B v v) {φ : E} (hφ : ‖φ‖ = 1) :
    RungStepCoupling (B φ φ) {bd | ∃ ψ : E, ‖ψ‖ = 1 ∧ ⟪φ, ψ⟫ = 0 ∧ bd = (B φ ψ, B ψ ψ)} := by
  rintro bd ⟨ψ, hψ, horth, rfl⟩
  have hc := coupling_le_of_margin B hB hμ hφ hψ horth
  have ha : μ ≤ B φ φ := by simpa [hφ] using hμ φ
  have hd : μ ≤ B ψ ψ := by simpa [hψ] using hμ ψ
  show (B φ ψ) ^ 2 ≤ B φ φ * B ψ ψ
  nlinarith

end GppWeilRungTwo

#print axioms GppWeilRungTwo.autocorr_neg
#print axioms GppWeilRungTwo.weilForm_rung_two
#print axioms GppWeilRungTwo.weilForm_rung_three
#print axioms GppWeilRungTwo.log_four_lt
#print axioms GppWeilRungTwo.rung_two_of_certificate
#print axioms GppWeilRungTwo.sq_le_mul_of_forall_nonneg
#print axioms GppWeilRungTwo.sq_le_of_margin
#print axioms GppWeilRungTwo.coupling_le_of_margin
#print axioms GppWeilRungTwo.rungStepCoupling_of_margin
