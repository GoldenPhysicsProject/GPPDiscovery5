/-
BTSigmaTranslationEvader.lean — GPPDiscovery5 (Heart of the Cosmos / Grok seat)
Task: Turbo, GPP-bridge PR #11 comment 6804647675 (7)(2): formalise the infinite-dimensional
one-parameter-group evader on L²(ℝ).

Content (all PROVED, standard axioms only, Mathlib-only imports):
* `transl W` : the translation isometry (V_W b)(x) = b(x + W) on Gp = L²(ℝ; ℂ), as a
  `LinearIsometry` (`Lp.compMeasurePreservingₗᵢ`); group law `transl_add`, `transl_zero`.
  Its generator −i d/dx is unbounded (not formalised; not needed).
* `inner_box_transl` : ⟪1_[0,1], V_W 1_[0,1]⟫ = max(0, 1 − |W|) exactly.
* `im_inner_evader` : with a = i·1_[0,1], b = 1_[0,1], Im⟪a, V_W b⟫ = −max(0, 1 − |W|).
  Via `blockForm_diag` (BTSigmaOneParam) this is the prime-block profile
  blockForm(a, V_W b, a, V_W b) = −2i·Im⟪a, V_W b⟫ (U_W = id, χ ≡ 1).
* `im_inner_evader_eq_zero_of_ge`, `im_inner_evader_zero` : vanishes for W ≥ 1, equals −1 at 0.
* `evader_not_analytic` : so W ↦ Im⟪a, V_W b⟫ is NOT real-analytic on ℝ — the hypothesis of
  `analytic_orbit_kill` (BTSigmaOneParam, 31ff892) genuinely fails for an unbounded-generator
  one-parameter group of isometries; the identity-theorem kill does not extend.

Scope / not claimed: the exact-profile realisation of an arbitrary C_c^∞ stratumSN via
Plancherel factorisation is paper level (not in this file). No zeros are used. S_b existence,
P-COV, P-UNIT, P-LOG, P-θ, P-LIM, BT-Σ, ARP and RH remain OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo.
-/
import Mathlib

open MeasureTheory Complex

noncomputable section

namespace BTSigmaTranslationEvader

/-- Gp = L²(ℝ; ℂ). -/
abbrev Gp := Lp ℂ 2 (volume : Measure ℝ)

/-- The translation isometry (V_W b)(x) = b(x + W). -/
def transl (W : ℝ) : Gp →ₗᵢ[ℂ] Gp :=
  Lp.compMeasurePreservingₗᵢ ℂ (fun x : ℝ => x + W) (measurePreserving_add_right volume W)

lemma transl_apply (W : ℝ) (f : Gp) :
    transl W f = Lp.compMeasurePreserving (fun x : ℝ => x + W)
      (measurePreserving_add_right volume W) f := rfl

/-- One-parameter group law V_{W+W'} = V_W ∘ V_{W'}. -/
theorem transl_add (W W' : ℝ) (f : Gp) : transl (W + W') f = transl W (transl W' f) := by
  rw [transl_apply, transl_apply, transl_apply]
  apply Lp.ext
  have h1 := Lp.coeFn_compMeasurePreserving f (measurePreserving_add_right volume (W + W'))
  have h2 := Lp.coeFn_compMeasurePreserving
    (Lp.compMeasurePreserving (fun x : ℝ => x + W') (measurePreserving_add_right volume W') f)
    (measurePreserving_add_right volume W)
  have h3 := Lp.coeFn_compMeasurePreserving f (measurePreserving_add_right volume W')
  have h3' := (measurePreserving_add_right volume W).quasiMeasurePreserving.ae_eq_comp h3
  filter_upwards [h1, h2, h3'] with x hx1 hx2 hx3
  rw [hx1, hx2, hx3]
  simp [Function.comp, add_assoc, add_comm W W']

theorem transl_zero (f : Gp) : transl 0 f = f := by
  rw [transl_apply]
  apply Lp.ext
  filter_upwards [Lp.coeFn_compMeasurePreserving f (measurePreserving_add_right volume 0)]
    with x hx
  rw [hx]; simp

/-- The unit box 1_[0,1] as a vector of Gp. -/
def box : Gp := indicatorConstLp 2 (measurableSet_Icc : MeasurableSet (Set.Icc (0:ℝ) 1))
  (by simp) (1 : ℂ)

lemma preimage_box (W : ℝ) :
    (fun x : ℝ => x + W) ⁻¹' Set.Icc (0:ℝ) 1 = Set.Icc (-W) (1 - W) := by
  ext x; simp only [Set.mem_preimage, Set.mem_Icc]
  constructor <;> rintro ⟨h1, h2⟩ <;> constructor <;> linarith

lemma volume_real_overlap (W : ℝ) :
    (volume : Measure ℝ).real (Set.Icc (0:ℝ) 1 ∩ Set.Icc (-W) (1 - W)) = max 0 (1 - |W|) := by
  rw [Set.Icc_inter_Icc, measureReal_def, Real.volume_Icc, ENNReal.toReal_ofReal']
  rcases le_total 0 W with h | h
  · rw [abs_of_nonneg h, max_eq_left (by linarith : -W ≤ 0), min_eq_right (by linarith),
      max_comm, sub_zero]
  · rw [abs_of_nonpos h, max_eq_right (by linarith : (0:ℝ) ≤ -W), min_eq_left (by linarith),
      max_comm, sub_neg_eq_add]

/-- The translation pairing ⟪1_[0,1], V_W 1_[0,1]⟫ = max(0, 1 − |W|) (the triangle profile). -/
theorem inner_box_transl (W : ℝ) :
    inner ℂ box (transl W box) = ((max 0 (1 - |W|) : ℝ) : ℂ) := by
  have hs : MeasurableSet ((fun x : ℝ => x + W) ⁻¹' Set.Icc (0:ℝ) 1) :=
    measurableSet_Icc.preimage (measurePreserving_add_right volume W).measurable
  rw [transl_apply, box, Lp.indicatorConstLp_compMeasurePreserving,
    MeasureTheory.L2.inner_indicatorConstLp_indicatorConstLp (hμs := _) (hμt := _), preimage_box,
    volume_real_overlap]
  simp

/-- The evader pair: a = i·1_[0,1], b = 1_[0,1]; Im⟪a, V_W b⟫ = −max(0, 1 − |W|). -/
theorem im_inner_evader (W : ℝ) :
    (inner ℂ (I • box) (transl W box)).im = -max 0 (1 - |W|) := by
  rw [inner_smul_left, inner_box_transl]
  simp

/-- The pairing vanishes past W = 1 … -/
theorem im_inner_evader_eq_zero_of_ge (W : ℝ) (hW : 1 ≤ W) :
    (inner ℂ (I • box) (transl W box)).im = 0 := by
  rw [im_inner_evader, abs_of_nonneg (by linarith), max_eq_left (by linarith)]; simp

/-- … without vanishing identically. -/
theorem im_inner_evader_zero : (inner ℂ (I • box) (transl 0 box)).im = -1 := by
  rw [im_inner_evader]; norm_num

/-- Hence the orbit pairing is NOT real-analytic on ℝ: the hypothesis of
`analytic_orbit_kill` (BTSigmaOneParam) fails for the translation group on L²(ℝ),
an unbounded-generator one-parameter group of isometries. -/
theorem evader_not_analytic :
    ¬ AnalyticOnNhd ℝ (fun W : ℝ => (inner ℂ (I • box) (transl W box)).im) Set.univ := by
  intro h
  have hev : (fun W : ℝ => (inner ℂ (I • box) (transl W box)).im) =ᶠ[nhds 2] 0 := by
    filter_upwards [Ioi_mem_nhds (by norm_num : (1:ℝ) < 2)] with W hW
    exact im_inner_evader_eq_zero_of_ge W (le_of_lt hW)
  have := h.eqOn_zero_of_preconnected_of_eventuallyEq_zero isPreconnected_univ
    (Set.mem_univ 2) hev (Set.mem_univ 0)
  simp only [Pi.zero_apply] at this
  rw [im_inner_evader_zero] at this
  norm_num at this

end BTSigmaTranslationEvader

#print axioms BTSigmaTranslationEvader.transl_add
#print axioms BTSigmaTranslationEvader.transl_zero
#print axioms BTSigmaTranslationEvader.inner_box_transl
#print axioms BTSigmaTranslationEvader.im_inner_evader
#print axioms BTSigmaTranslationEvader.im_inner_evader_eq_zero_of_ge
#print axioms BTSigmaTranslationEvader.im_inner_evader_zero
#print axioms BTSigmaTranslationEvader.evader_not_analytic
