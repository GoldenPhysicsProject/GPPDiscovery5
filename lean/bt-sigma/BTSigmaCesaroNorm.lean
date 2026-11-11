/-
BTSigmaCesaroNorm.lean — Turbo 6803068030 (Daniel's ruling: statehood norm = Cesàro mean).

In the log-norm coordinate u = log |x| on the norm direction of the idele class group, a
character of weight s = 1/2 + δ + iγ has norm density |e^{(s-1/2)u}|² = e^{2δu}.  The symmetric
Cesàro mean over [-T, T] is
    cesaroMean δ T = (1/(2T)) ∫_{-T}^{T} e^{2δu} du.
We prove: it equals sinh(2δT)/(2δT) for δ ≠ 0, equals 1 for δ = 0 (every T > 0), is even in δ
(the shadow involution u ↦ -u swaps the two ends), and tends to +∞ as T → ∞ for every δ ≠ 0.
So "Cesàro-finite" (bounded in T) ⇔ δ = 0, with limit norm 1 on the line.

Labels: all theorems PROVED (standard axioms).  Scope: scalar density only; nothing here
constructs the boundary embedding or S_b (OPEN), and RH is OPEN.
Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok). Task: Turbo.
-/
import Mathlib

open Real MeasureTheory Filter Topology intervalIntegral

namespace BTSigmaCesaro

/-- Symmetric Cesàro mean of the norm density e^{2δu} over [-T, T]. -/
noncomputable def cesaroMean (δ T : ℝ) : ℝ :=
  (1 / (2 * T)) * ∫ u in (-T)..T, Real.exp (2 * δ * u)

/-- Off the line: closed form sinh(2δT)/(2δT). -/
theorem cesaroMean_eq_sinh {δ T : ℝ} (hδ : δ ≠ 0) (hT : T ≠ 0) :
    cesaroMean δ T = Real.sinh (2 * δ * T) / (2 * δ * T) := by
  unfold cesaroMean
  have h2 : (2 * δ) ≠ 0 := by positivity
  have hint : ∫ u in (-T)..T, Real.exp (2 * δ * u)
      = (Real.exp (2 * δ * T) - Real.exp (2 * δ * (-T))) / (2 * δ) := by
    rw [intervalIntegral.integral_comp_mul_left (fun x => Real.exp x) h2, integral_exp, smul_eq_mul]
    field_simp
  rw [hint, Real.sinh_eq]
  have : 2 * δ * (-T) = -(2 * δ * T) := by ring
  rw [this]
  field_simp

/-- On the line: exactly 1 for every T > 0. -/
theorem cesaroMean_zero {T : ℝ} (hT : 0 < T) : cesaroMean 0 T = 1 := by
  unfold cesaroMean
  simp
  field_simp
  ring

/-- Shadow involution: the symmetric mean is even in δ. -/
theorem cesaroMean_neg (δ T : ℝ) : cesaroMean (-δ) T = cesaroMean δ T := by
  unfold cesaroMean
  congr 1
  have := intervalIntegral.integral_comp_neg (a := -T) (b := T)
    (fun u => Real.exp (2 * δ * u))
  simp only [neg_neg] at this
  rw [← this]
  congr 1; funext u; ring_nf

/-- Elementary lower bound: sinh x / x ≥ (1 + x/2)/2 for x > 0. -/
theorem sinh_div_ge {x : ℝ} (hx : 0 < x) : (1 + x / 2) / 2 ≤ Real.sinh x / x := by
  have h1 := Real.quadratic_le_exp_of_nonneg hx.le
  have h2 : Real.exp (-x) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
  rw [le_div_iff₀ hx, Real.sinh_eq]
  nlinarith

/-- Off the line (δ > 0): linear lower bound in T. -/
theorem cesaroMean_ge_of_pos {δ T : ℝ} (hδ : 0 < δ) (hT : 0 < T) :
    (1 + δ * T) / 2 ≤ cesaroMean δ T := by
  rw [cesaroMean_eq_sinh hδ.ne' hT.ne']
  have := sinh_div_ge (x := 2 * δ * T) (by positivity)
  have e : (1 + 2 * δ * T / 2) / 2 = (1 + δ * T) / 2 := by ring
  linarith

/-- Off the line: the Cesàro norm diverges for every δ ≠ 0. -/
theorem cesaroMean_tendsto_atTop {δ : ℝ} (hδ : δ ≠ 0) :
    Tendsto (fun T => cesaroMean δ T) atTop atTop := by
  have hpos : ∀ {d : ℝ}, 0 < d → Tendsto (fun T => cesaroMean d T) atTop atTop := by
    intro d hd
    have hl : Tendsto (fun T : ℝ => (1 + d * T) / 2) atTop atTop := by
      apply Tendsto.atTop_div_const (by norm_num)
      apply tendsto_atTop_add_const_left
      exact Tendsto.const_mul_atTop hd tendsto_id
    refine tendsto_atTop_mono' atTop ?_ hl
    filter_upwards [eventually_gt_atTop 0] with T hT
    exact cesaroMean_ge_of_pos hd hT
  rcases lt_or_gt_of_ne hδ with h | h
  · have := hpos (neg_pos.mpr h)
    simpa [cesaroMean_neg] using this
  · exact hpos h

/-- Statehood: a weight is a generalized state iff its Cesàro norm stays bounded. -/
def CesaroFinite (δ : ℝ) : Prop := ∃ C : ℝ, ∀ T : ℝ, 0 < T → cesaroMean δ T ≤ C

/-- Cesàro-finite exactly on the principal series line δ = 0. -/
theorem cesaroFinite_iff (δ : ℝ) : CesaroFinite δ ↔ δ = 0 := by
  constructor
  · rintro ⟨C, hC⟩
    by_contra hδ
    have ht := (cesaroMean_tendsto_atTop hδ).eventually (eventually_gt_atTop C)
    obtain ⟨T, hT1, hT2⟩ := (ht.and (eventually_gt_atTop 0)).exists
    exact absurd (hC T hT2) (not_le.mpr hT1)
  · rintro rfl
    exact ⟨1, fun T hT => (cesaroMean_zero hT).le⟩

/-- Named failure (Turbo 6803068030): any norm assignment that is finite at some δ ≠ 0
is not the Cesàro norm there. -/
theorem finite_offline_not_cesaro {δ : ℝ} (hδ : δ ≠ 0) : ¬ CesaroFinite δ := by
  rw [cesaroFinite_iff]; exact hδ

end BTSigmaCesaro

#print axioms BTSigmaCesaro.cesaroMean_eq_sinh
#print axioms BTSigmaCesaro.cesaroMean_zero
#print axioms BTSigmaCesaro.cesaroMean_neg
#print axioms BTSigmaCesaro.cesaroMean_tendsto_atTop
#print axioms BTSigmaCesaro.cesaroFinite_iff
#print axioms BTSigmaCesaro.finite_offline_not_cesaro
