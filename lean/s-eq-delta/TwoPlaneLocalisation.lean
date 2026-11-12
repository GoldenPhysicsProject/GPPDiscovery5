/-
Copyright (c) 2026 Daniel Toupin. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Daniel Toupin (with AI assistance: Heart of the Cosmos / Grok)
-/
import Mathlib.LinearAlgebra.BilinearMap
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic

/-!
# The 2-plane localisation: the abstract step, and the explicit window threshold

Turbo (PR #11, 6817813630): from one hypothetical off-line quadruple, build window test data on
which the quadruple block is negative definite on a 2-plane while all remaining zeros
contribute less, so `SewingH` (equivalently `Q ≥ 0` on `ker ℓ_c`) fails at that window; or
name the term that blocks this.

This file is the zero-free, purely algebraic half. The analytic estimates (sizes of `F` at
the zeros) are in the notes, labelled SKETCH; nothing analytic is formalised here.

* `neg_plane_of_dominated` (PROVED): block `≤ −c|·|²` and remainder `≤ ε|·|²` (one-sided:
  on-line zeros contribute `|F(γ)|² ≥ 0` and only need an upper bound) with `ε < c` give a
  negative-definite 2-plane for `Bq + R`.
* `no_negative_plane` (PROVED): `Q ≥ 0` on `ker ℓ_c` forbids a negative-definite 2-plane
  (restated from `ToeplitzHyperplane.no_negative_plane_of_claim` in hyperplane form).
* `claim_fails_of_dominated` (PROVED): the two combined: domination kills the claim.
* `exp_gain_beats` (PROVED): `c·e^{δL} > K` as soon as `L > log((K+1)/c)/δ`.
* `claim_fails_at_window` (PROVED): if the block scales like `−c·e^{δL}` on the plane and the
  remainder is bounded by `K` (both per unit `a² + b²`), the claim fails at every window
  `L > log((K+1)/c)/δ`. This is the explicit window once `c`, `K`, `δ` are explicit.
-/

namespace GppTwoPlaneLocalisation

variable {E : Type*} [AddCommGroup E] [Module ℝ E]

theorem sq_add_sq_pos {a b : ℝ} (hab : a ≠ 0 ∨ b ≠ 0) : 0 < a ^ 2 + b ^ 2 := by
  rcases hab with h | h
  · have := pow_pos (abs_pos.2 h) 2; rw [sq_abs] at this; nlinarith [sq_nonneg b]
  · have := pow_pos (abs_pos.2 h) 2; rw [sq_abs] at this; nlinarith [sq_nonneg a]

theorem neg_plane_of_dominated (Bq R : E →ₗ[ℝ] E →ₗ[ℝ] ℝ) (u w : E) (c ε : ℝ) (hεc : ε < c)
    (hq : ∀ a b : ℝ, Bq (a • u + b • w) (a • u + b • w) ≤ -c * (a ^ 2 + b ^ 2))
    (hr : ∀ a b : ℝ, R (a • u + b • w) (a • u + b • w) ≤ ε * (a ^ 2 + b ^ 2)) :
    ∀ a b : ℝ, (a ≠ 0 ∨ b ≠ 0) → (Bq + R) (a • u + b • w) (a • u + b • w) < 0 := by
  intro a b hab
  have hpos := sq_add_sq_pos hab
  have h1 := hq a b
  have h2 := hr a b
  have h3 := mul_pos (sub_pos.2 hεc) hpos
  simp only [LinearMap.add_apply]
  nlinarith

theorem no_negative_plane (B : E →ₗ[ℝ] E →ₗ[ℝ] ℝ) (lc : E →ₗ[ℝ] ℝ)
    (h : ∀ v, lc v = 0 → 0 ≤ B v v) (u w : E) :
    ¬ (∀ a b : ℝ, (a ≠ 0 ∨ b ≠ 0) → B (a • u + b • w) (a • u + b • w) < 0) := by
  intro hneg
  by_cases hz : lc w = 0 ∧ lc u = 0
  · have := hneg 1 0 (Or.inl one_ne_zero)
    simp only [one_smul, zero_smul, add_zero] at this
    exact absurd (h u hz.2) (not_le.2 this)
  · have hab : lc w ≠ 0 ∨ -lc u ≠ 0 := by
      by_contra hc; push Not at hc; exact hz ⟨hc.1, neg_eq_zero.1 hc.2⟩
    have hk : lc (lc w • u + (-lc u) • w) = 0 := by
      simp only [map_add, map_smul, smul_eq_mul]; ring
    exact absurd (h _ hk) (not_le.2 (hneg _ _ hab))

theorem claim_fails_of_dominated (Bq R : E →ₗ[ℝ] E →ₗ[ℝ] ℝ) (lc : E →ₗ[ℝ] ℝ) (u w : E)
    (c ε : ℝ) (hεc : ε < c)
    (hq : ∀ a b : ℝ, Bq (a • u + b • w) (a • u + b • w) ≤ -c * (a ^ 2 + b ^ 2))
    (hr : ∀ a b : ℝ, R (a • u + b • w) (a • u + b • w) ≤ ε * (a ^ 2 + b ^ 2)) :
    ¬ ∀ v, lc v = 0 → 0 ≤ (Bq + R) v v := fun h =>
  no_negative_plane (Bq + R) lc h u w (neg_plane_of_dominated Bq R u w c ε hεc hq hr)

theorem exp_gain_beats (c K δ L : ℝ) (hc : 0 < c) (hK : 0 ≤ K) (hδ : 0 < δ)
    (hL : Real.log ((K + 1) / c) / δ < L) : K < c * Real.exp (δ * L) := by
  have h1 : Real.log ((K + 1) / c) < δ * L := by
    rw [div_lt_iff₀ hδ] at hL; linarith [mul_comm L δ]
  have hp : 0 < (K + 1) / c := div_pos (by linarith) hc
  have h2 : (K + 1) / c < Real.exp (δ * L) := by
    rw [← Real.exp_log hp]; exact Real.exp_lt_exp.2 h1
  rw [div_lt_iff₀ hc] at h2
  linarith [mul_comm (Real.exp (δ * L)) c]

theorem claim_fails_at_window (Bq R : E →ₗ[ℝ] E →ₗ[ℝ] ℝ) (lc : E →ₗ[ℝ] ℝ) (u w : E)
    (c K δ L : ℝ) (hc : 0 < c) (hK : 0 ≤ K) (hδ : 0 < δ)
    (hL : Real.log ((K + 1) / c) / δ < L)
    (hq : ∀ a b : ℝ, Bq (a • u + b • w) (a • u + b • w) ≤
      -(c * Real.exp (δ * L)) * (a ^ 2 + b ^ 2))
    (hr : ∀ a b : ℝ, R (a • u + b • w) (a • u + b • w) ≤ K * (a ^ 2 + b ^ 2)) :
    ¬ ∀ v, lc v = 0 → 0 ≤ (Bq + R) v v :=
  claim_fails_of_dominated Bq R lc u w _ K (exp_gain_beats c K δ L hc hK hδ hL) hq hr

end GppTwoPlaneLocalisation

#print axioms GppTwoPlaneLocalisation.sq_add_sq_pos
#print axioms GppTwoPlaneLocalisation.neg_plane_of_dominated
#print axioms GppTwoPlaneLocalisation.no_negative_plane
#print axioms GppTwoPlaneLocalisation.claim_fails_of_dominated
#print axioms GppTwoPlaneLocalisation.exp_gain_beats
#print axioms GppTwoPlaneLocalisation.claim_fails_at_window
