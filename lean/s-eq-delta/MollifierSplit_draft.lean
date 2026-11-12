/-
Copyright (c) 2026 GoldenPhysicsProject. Released under the project license.
Heart of the Cosmos (Grok, GPPDiscovery5 seat). Turbo PR #11 6827581325, "next exact instruction".

# The shared mollification root of `Regularity` and `ArchMatch`, split

`g_d = ½[triG(· − d) + triG(· + d)]` is `C²_c`, not `C^∞_c`, so the in-tree explicit formula
(`explicit_formula_even`) and Archimedean identification (`archimedean_identification`) do not
apply to it directly. Mollify: `φ_n = (bump n).normed` (Mathlib `ContDiffBump` at `0`,
`rOut = 1/(n+2)`), `g_{d,n} = φ_n ⋆ g_d`.

PROVED from Mathlib / in tree: `moll_isTest` (`g_{d,n}` is `C^∞_c`, `contDiff_convolution_left`),
`moll_even`, `moll_eq_zero` (support in `[−|d|−3, |d|+3]`), `abs_moll_le`, `regSumm` (the
summability half of `Regularity`, unconditionally), `primeSide_tendsto` (finite sum),
`poleTerm_tendsto` (dominated convergence).

Named Props, stated exactly (NOT proved):
* (a) `MollC1`: `g_{d,n} → g_d` and `g_{d,n}' → g_d'` uniformly.
* (b) `StripFactor`: on `0 ≤ Re s ≤ 1`, `ĝ_{d,n}(s) = ĝ_d(s) φ̂_n(s)`, `|φ̂_n(s)| ≤ e^{1/2}`,
  and `φ̂_n(s) → 1`.
* (c) `RealPlaceDCT`: dominated convergence for the real-place integral on `(0,∞)` under C¹
  convergence with uniform support.

PROVED conditional assemblies: `regularity_of_moll : (a) → (b) → (c) → Regularity`,
`archMatch_of_moll : (a) → (b) → (c) → ArchMatch` (via `archimedean_identification` for each
`g_{d,n}`, the limit `archId_gd`, Fourier inversion at 0 for `g_d` by the same limit, the identity
`ĝ_d(1/2 + ir) = Ψ(r)² cos(rd)`, `ψ(z̄) = conj ψ(z)` and the even fold `integral_comp_abs`),
`hefIdent_of_moll`, `hPsiMollChain : HPsi → (a) → (b) → (c) → KreinDecomp → IdentityStep → RH`.
H_ψ: OPEN. ARP: OPEN. RH: OPEN.
-/
import Mathlib.Analysis.Normed.Group.Tannery
import GppVerify.RiemannHypothesis.HEFIdentSplit
import GppVerify.RiemannHypothesis.ArchimedeanIdentification
import GppVerify.RiemannHypothesis.ShiftContour
import GppVerify.RiemannHypothesis.DigammaQuarterSech
import GppVerify.RiemannHypothesis.DigammaGaussSeries

namespace GppMollifierSplit

open scoped Convolution
open MeasureTheory Filter Topology GppKreinHPsi GppKreinHPsiChain GppKreinPoleSplit
  GppKreinHEFSplit GppHEFIdentSplit GppKreinRepSplit

/-! ## 0. The concrete mollifier `φ_n` and `g_{d,n} = φ_n ⋆ g_d` -/

/-- Bump at `0` with `rIn = 1/(2(n+2))`, `rOut = 1/(n+2)` (so `ε_n = 1/(n+2) ≤ 1/2`). -/
noncomputable def bump (n : ℕ) : ContDiffBump (0 : ℝ) where
  rIn := 1 / (2 * ((n : ℝ) + 2))
  rOut := 1 / ((n : ℝ) + 2)
  rIn_pos := by positivity
  rIn_lt_rOut := by
    have h : (0 : ℝ) < (n : ℝ) + 2 := by positivity
    apply one_div_lt_one_div_of_lt h
    linarith

/-- The mollifier `g_{d,n} = φ_n ⋆ g_d`, `φ_n = (bump n).normed` (even, `≥ 0`, `∫ = 1`, smooth). -/
noncomputable def moll (d : ℝ) (n : ℕ) : ℝ → ℝ :=
  (bump n).normed volume ⋆[ContinuousLinearMap.lsmul ℝ ℝ, volume] gd d

theorem moll_apply (d : ℝ) (n : ℕ) (x : ℝ) :
    moll d n x = ∫ t, (bump n).normed volume t * gd d (x - t) := by
  rw [moll, convolution_lsmul]; rfl

theorem bump_eq_zero {n : ℕ} {t : ℝ} (ht : 1 / 2 ≤ |t|) : (bump n).normed volume t = 0 := by
  by_contra h
  have hm : t ∈ Function.support ((bump n).normed volume) := h
  rw [(bump n).support_normed_eq, Metric.mem_ball, Real.dist_eq, sub_zero] at hm
  have hr : (bump n).rOut ≤ 1 / 2 := by
    show 1 / ((n : ℝ) + 2) ≤ 1 / 2
    exact one_div_le_one_div_of_le (by norm_num) (by linarith [Nat.cast_nonneg (α := ℝ) n])
  linarith

theorem gd_eq_zero {d u : ℝ} (h : |d| + 2 < |u|) : gd d u = 0 := by
  unfold gd
  have h1 : |u| - |d| ≤ |u - d| := abs_sub_abs_le_abs_sub u d
  have h2 : |u| - |d| ≤ |u + d| := by
    have := abs_sub_abs_le_abs_sub u (-d); rwa [abs_neg, sub_neg_eq_add] at this
  rw [triG_eq_zero (by linarith), triG_eq_zero (by linarith)]; ring

theorem gd_measurable (d : ℝ) : Measurable (gd d) := by
  unfold gd
  exact ((triG_measurable.comp (measurable_id.sub_const d)).add
    (triG_measurable.comp (measurable_id.add_const d))).div_const 2

theorem abs_gd_le (d u : ℝ) : |gd d u| ≤ 1 := by
  unfold gd
  have h1 := abs_le.1 (abs_triG_le (u - d))
  have h2 := abs_le.1 (abs_triG_le (u + d))
  rw [abs_le]; constructor <;> linarith [h1.1, h1.2, h2.1, h2.2]

theorem gd_integrable (d : ℝ) : Integrable (gd d) := by
  have hb : Integrable ((Set.Icc (-(|d| + 2)) (|d| + 2)).indicator fun _ => (1 : ℝ)) :=
    (integrable_indicator_iff measurableSet_Icc).2 (integrableOn_const (by simp))
  refine hb.mono' (gd_measurable d).aestronglyMeasurable (Eventually.of_forall fun u => ?_)
  by_cases hu : u ∈ Set.Icc (-(|d| + 2)) (|d| + 2)
  · rw [Set.indicator_of_mem hu, Real.norm_eq_abs]; exact abs_gd_le d u
  · rw [Set.indicator_of_notMem hu, Real.norm_eq_abs, gd_eq_zero, abs_zero]
    rw [Set.mem_Icc, not_and_or, not_le, not_le] at hu
    have hd := abs_nonneg d
    rcases hu with hu | hu
    · rw [abs_of_neg (show u < 0 by linarith)]; linarith
    · rw [abs_of_pos (show 0 < u by linarith)]; linarith

theorem moll_eq_zero {d u : ℝ} (n : ℕ) (h : |d| + 3 < |u|) : moll d n u = 0 := by
  rw [moll_apply]
  have hz : ∀ t, (bump n).normed volume t * gd d (u - t) = 0 := by
    intro t
    by_cases ht : 1 / 2 ≤ |t|
    · rw [bump_eq_zero ht, zero_mul]
    · push Not at ht
      have : |u| - |t| ≤ |u - t| := abs_sub_abs_le_abs_sub u t
      rw [gd_eq_zero (by linarith), mul_zero]
  simp only [hz, integral_zero]

/-- `g_{d,n}` is a test function (`C^∞_c`): Mathlib `contDiff_convolution_left`. PROVED. -/
theorem moll_isTest (d : ℝ) (n : ℕ) : GppMeyerSpectralWeil.IsTestFunction (moll d n) := by
  refine ⟨?_, ?_⟩
  · exact (bump n).hasCompactSupport_normed.contDiff_convolution_left (ContinuousLinearMap.lsmul ℝ ℝ) (bump n).contDiff_normed
      (gd_integrable d).locallyIntegrable
  · refine HasCompactSupport.intro (isCompact_Icc (a := -(|d| + 3)) (b := |d| + 3)) fun u hu => ?_
    apply moll_eq_zero
    rw [Set.mem_Icc, not_and_or, not_le, not_le] at hu
    have hd := abs_nonneg d
    rcases hu with hu | hu
    · rw [abs_of_neg (show u < 0 by linarith)]; linarith
    · rw [abs_of_pos (show 0 < u by linarith)]; linarith

/-- `g_{d,n}` is even. PROVED. -/
theorem moll_even (d : ℝ) (n : ℕ) (u : ℝ) : moll d n (-u) = moll d n u := by
  rw [moll_apply, moll_apply,
    ← integral_neg_eq_self (fun t => (bump n).normed volume t * gd d (-u - t)) volume]
  congr 1; funext t
  rw [(bump n).normed_neg, show -u - -t = -(u - t) by ring, gd_even]

/-- `|g_{d,n}| ≤ 1`. PROVED. -/
theorem abs_moll_le (d : ℝ) (n : ℕ) (u : ℝ) : |moll d n u| ≤ 1 := by
  rw [moll_apply, ← Real.norm_eq_abs]
  calc ‖∫ t, (bump n).normed volume t * gd d (u - t)‖ ≤ ∫ t, (bump n).normed volume t :=
        norm_integral_le_of_norm_le (bump n).integrable_normed (Eventually.of_forall fun t => by
          rw [norm_mul, Real.norm_eq_abs, Real.norm_eq_abs,
            abs_of_nonneg ((bump n).nonneg_normed t)]
          exact mul_le_of_le_one_right ((bump n).nonneg_normed t) (abs_gd_le d (u - t)))
    _ = 1 := (bump n).integral_normed

/-! ## 1. The three missing lemmas, stated exactly as named Props -/

/-- **(a) `MollC1`** (named premise): C¹ convergence of the mollifiers of the `C¹_c` function
`g_d`: `g_{d,n} → g_d` and `g_{d,n}' → g_d'` uniformly on `ℝ` (all supports lie in the fixed
compact `[−|d|−3, |d|+3]`, so this is C¹ convergence on the needed compacta). -/
def MollC1 : Prop := ∀ d : ℝ,
  TendstoUniformly (fun n => moll d n) (gd d) atTop ∧
    TendstoUniformly (fun n => deriv (moll d n)) (deriv (gd d)) atTop

/-- `φ̂_n(s) = ∫ φ_n(u) e^{(s − 1/2)u} du`. -/
noncomputable def mhat (n : ℕ) (s : ℂ) : ℂ :=
  GppMeyerSpectralWeil.spectralTransform ((bump n).normed volume) s

/-- **(b) `StripFactor`** (named premise): on the closed strip `0 ≤ Re s ≤ 1`, the convolution
/ Fubini identity `ĝ_{d,n}(s) = ĝ_d(s) φ̂_n(s)`, the strip bound `|φ̂_n(s)| ≤ e^{1/2}`
(true bound `≤ e^{ε_n/2}`), and the pointwise limit `φ̂_n(s) → 1`. -/
def StripFactor : Prop :=
  (∀ (d : ℝ) (n : ℕ) (s : ℂ), 0 ≤ s.re → s.re ≤ 1 →
      GppMeyerSpectralWeil.spectralTransform (moll d n) s =
        GppMeyerSpectralWeil.spectralTransform (gd d) s * mhat n s) ∧
    (∀ (n : ℕ) (s : ℂ), 0 ≤ s.re → s.re ≤ 1 → ‖mhat n s‖ ≤ Real.exp (1 / 2)) ∧
    (∀ s : ℂ, 0 ≤ s.re → s.re ≤ 1 → Tendsto (fun n => mhat n s) atTop (𝓝 1))

/-- **(c) `RealPlaceDCT`** (named premise): dominated convergence for the real-place integral on
`(0,∞)` under C¹ convergence with uniform support. -/
def RealPlaceDCT : Prop := ∀ (G : ℕ → ℝ → ℝ) (g : ℝ → ℝ) (R : ℝ),
  (∀ n, ContDiff ℝ 1 (G n)) → (∀ n u, R < |u| → G n u = 0) → (∀ u, R < |u| → g u = 0) →
  TendstoUniformly G g atTop → TendstoUniformly (fun n => deriv (G n)) (deriv g) atTop →
  Tendsto (fun n => ∫ u in Set.Ioi (0 : ℝ), GppArchimedeanFloor.archimedeanIntegrand (G n) u)
    atTop (𝓝 (∫ u in Set.Ioi (0 : ℝ), GppArchimedeanFloor.archimedeanIntegrand g u))

/-! ## 2. Unconditional pieces -/

/-- **Summability half of `Regularity`, PROVED unconditionally** (`cosToExp` + `transform`, index
transport `toEquiv`). -/
theorem regSumm (d : ℝ) : Summable (fun ρ : GppWeilCriterion.nontrivialZeros =>
    (GppXiZeroCounting.zeroMult ρ : ℂ) * GppMeyerSpectralWeil.spectralTransform (gd d) ρ) := by
  rw [← toEquiv.summable_iff]
  have h := (cosToExp d).1
  unfold ZeroSummable at h
  refine h.congr fun γ => ?_
  show (mult (γ : ℂ) : ℂ) * (efCoeff γ * Complex.cos (γ * d)) =
    (GppXiZeroCounting.zeroMult (rho (γ : ℂ)) : ℂ) *
      GppMeyerSpectralWeil.spectralTransform (gd d) (rho γ)
  rw [transform d γ γ.2]; rfl

/-- Prime side: finite sum, pointwise convergence. PROVED. -/
theorem primeSide_tendsto {G : ℕ → ℝ → ℝ} {g : ℝ → ℝ} {R : ℝ}
    (hG : ∀ n u, R < |u| → G n u = 0) (hg : ∀ u, R < |u| → g u = 0)
    (hpt : ∀ u, Tendsto (fun n => G n u) atTop (𝓝 (g u))) :
    Tendsto (fun n => GppWeilLadder.primeSide (G n)) atTop (𝓝 (GppWeilLadder.primeSide g)) := by
  set N : ℕ := ⌈Real.exp R⌉₊ + 1 with hN
  have hbig : ∀ m : ℕ, m ∉ Finset.range N → R < |Real.log m| := by
    intro m hm
    rw [Finset.mem_range, not_lt] at hm
    have h1 : Real.exp R < m := by
      have hc := Nat.le_ceil (Real.exp R)
      have : ((⌈Real.exp R⌉₊ + 1 : ℕ) : ℝ) ≤ m := by exact_mod_cast hm
      push_cast at this; linarith
    have hm0 : (0 : ℝ) < m := lt_trans (Real.exp_pos R) h1
    have hm1 : (1 : ℝ) ≤ m := by
      have : 1 ≤ m := le_trans (Nat.le_add_left 1 _) hm
      exact_mod_cast this
    rw [abs_of_nonneg (Real.log_nonneg hm1)]
    exact (Real.lt_log_iff_exp_lt hm0).2 h1
  have hterm : ∀ f : ℝ → ℝ, (∀ u, R < |u| → f u = 0) → ∀ m ∉ Finset.range N,
      ArithmeticFunction.vonMangoldt m / Real.sqrt m * (f (Real.log m) + f (-(Real.log m))) = 0 := by
    intro f hf m hm
    have h := hbig m hm
    rw [hf _ h, hf _ (by rwa [abs_neg]), add_zero, mul_zero]
  unfold GppWeilLadder.primeSide
  have e : (fun n => ∑' m : ℕ, ArithmeticFunction.vonMangoldt m / Real.sqrt m *
      (G n (Real.log m) + G n (-(Real.log m)))) = fun n => ∑ m ∈ Finset.range N,
      ArithmeticFunction.vonMangoldt m / Real.sqrt m * (G n (Real.log m) + G n (-(Real.log m))) :=
    funext fun n => tsum_eq_sum (hterm _ (hG n))
  rw [e, tsum_eq_sum (hterm g hg)]
  exact tendsto_finsetSum _ fun m _ => ((hpt _).add (hpt _)).const_mul _

/-- Pole term: dominated convergence with the bound `1_{[−R,R]} · 2cosh(u/2)`. PROVED. -/
theorem poleTerm_tendsto (d : ℝ) (hU : TendstoUniformly (fun n => moll d n) (gd d) atTop) :
    Tendsto (fun n => GppArithmeticReflection.poleTerm (moll d n)) atTop
      (𝓝 (GppArithmeticReflection.poleTerm (gd d))) := by
  unfold GppArithmeticReflection.poleTerm
  set R : ℝ := |d| + 3
  refine tendsto_integral_of_dominated_convergence
    ((Set.Icc (-R) R).indicator fun u => 2 * Real.cosh (u / 2)) (fun n => ?_) ?_ (fun n => ?_)
    (Eventually.of_forall fun u => ?_)
  · exact ((moll_isTest d n).1.continuous.mul (by fun_prop)).aestronglyMeasurable
  · exact (integrable_indicator_iff measurableSet_Icc).2
      (Continuous.integrableOn_Icc (by fun_prop))
  · refine Eventually.of_forall fun u => ?_
    by_cases hu : u ∈ Set.Icc (-R) R
    · rw [Set.indicator_of_mem hu, norm_mul, Real.norm_eq_abs, Real.norm_eq_abs,
        abs_of_pos (by positivity : (0 : ℝ) < 2 * Real.cosh (u / 2))]
      exact mul_le_of_le_one_left (by positivity) (abs_moll_le d n u)
    · rw [Set.indicator_of_notMem hu, moll_eq_zero n, zero_mul, norm_zero]
      rw [Set.mem_Icc, not_and_or, not_le, not_le] at hu
      have hd := abs_nonneg d
      rcases hu with hu | hu
      · rw [abs_of_neg (show u < 0 by linarith)]; linarith
      · rw [abs_of_pos (show 0 < u by linarith)]; linarith
  · exact (hU.tendsto_at u).mul_const _

/-- Real-place term: `(a)` + `(c)`. PROVED (conditional on the named Props). -/
theorem archLimit (ha : MollC1) (hc : RealPlaceDCT) (d : ℝ) :
    Tendsto (fun n => GppArchimedeanFloor.archimedeanTerm (moll d n)) atTop
      (𝓝 (GppArchimedeanFloor.archimedeanTerm (gd d))) := by
  obtain ⟨hU, hD⟩ := ha d
  have hsuppM : ∀ n u, |d| + 3 < |u| → moll d n u = 0 := fun n u h => moll_eq_zero n h
  have hsuppG : ∀ u, |d| + 3 < |u| → gd d u = 0 := fun u h => gd_eq_zero (by linarith)
  unfold GppArchimedeanFloor.archimedeanTerm
  exact ((hU.tendsto_at 0).const_mul _).add ((hc (moll d) (gd d) (|d| + 3)
    (fun n => (moll_isTest d n).1.of_le (by exact_mod_cast le_top)) hsuppM hsuppG hU hD).const_mul 2)

/-! ## 3. Assembly: `(a) + (b) + (c) → Regularity` (PROVED) -/

theorem regularity_of_moll (ha : MollC1) (hb : StripFactor) (hc : RealPlaceDCT) :
    Regularity := by
  intro d
  refine ⟨regSumm d, ?_⟩
  obtain ⟨hU, hD⟩ := ha d
  obtain ⟨hfac, hbd, hlim⟩ := hb
  have hEF : ∀ n : ℕ, ∑' ρ : GppWeilCriterion.nontrivialZeros,
      (GppXiZeroCounting.zeroMult ρ : ℂ) * GppMeyerSpectralWeil.spectralTransform (moll d n) ρ =
        (GppArithmeticReflection.weilForm (moll d n) : ℂ) := fun n =>
    (GppExplicitFormulaAssembly.explicit_formula_even (moll_isTest d n) (moll_even d n)).2
  have hstrip : ∀ ρ : GppWeilCriterion.nontrivialZeros, 0 ≤ (ρ : ℂ).re ∧ (ρ : ℂ).re ≤ 1 :=
    fun ρ => ⟨ρ.2.2.1.le, ρ.2.2.2.le⟩
  have hZ : Tendsto (fun n => ∑' ρ : GppWeilCriterion.nontrivialZeros,
      (GppXiZeroCounting.zeroMult ρ : ℂ) * GppMeyerSpectralWeil.spectralTransform (moll d n) ρ)
      atTop (𝓝 (∑' ρ : GppWeilCriterion.nontrivialZeros,
      (GppXiZeroCounting.zeroMult ρ : ℂ) * GppMeyerSpectralWeil.spectralTransform (gd d) ρ)) := by
    refine tendsto_tsum_of_dominated_convergence (bound := fun ρ => Real.exp (1 / 2) *
      ‖(GppXiZeroCounting.zeroMult ρ : ℂ) * GppMeyerSpectralWeil.spectralTransform (gd d) ρ‖)
      ((summable_norm_iff.2 (regSumm d)).mul_left _) (fun ρ => ?_)
      (Eventually.of_forall fun n ρ => ?_)
    · have e : ∀ n : ℕ, (GppXiZeroCounting.zeroMult ρ : ℂ) *
          GppMeyerSpectralWeil.spectralTransform (moll d n) ρ =
          ((GppXiZeroCounting.zeroMult ρ : ℂ) *
            GppMeyerSpectralWeil.spectralTransform (gd d) ρ) * mhat n ρ := fun n => by
        rw [hfac d n ρ (hstrip ρ).1 (hstrip ρ).2]; ring
      simp_rw [e]
      simpa using (hlim ρ (hstrip ρ).1 (hstrip ρ).2).const_mul
        ((GppXiZeroCounting.zeroMult ρ : ℂ) * GppMeyerSpectralWeil.spectralTransform (gd d) ρ)
    · rw [hfac d n ρ (hstrip ρ).1 (hstrip ρ).2, ← mul_assoc, norm_mul, mul_comm (Real.exp _)]
      exact mul_le_mul_of_nonneg_left (hbd n ρ (hstrip ρ).1 (hstrip ρ).2) (norm_nonneg _)
  set R : ℝ := |d| + 3 with hR
  have hsuppM : ∀ n u, R < |u| → moll d n u = 0 := fun n u h => moll_eq_zero n h
  have hsuppG : ∀ u, R < |u| → gd d u = 0 := fun u h => gd_eq_zero (by linarith)
  have hpt : ∀ u, Tendsto (fun n => moll d n u) atTop (𝓝 (gd d u)) := fun u => hU.tendsto_at u
  have hA := archLimit ha hc d
  have hW : Tendsto (fun n => (GppArithmeticReflection.weilForm (moll d n) : ℂ)) atTop
      (𝓝 (GppArithmeticReflection.weilForm (gd d) : ℂ)) := by
    have h : Tendsto (fun n => GppArithmeticReflection.weilForm (moll d n)) atTop
        (𝓝 (GppArithmeticReflection.weilForm (gd d))) := by
      unfold GppArithmeticReflection.weilForm
      exact ((poleTerm_tendsto d hU).add hA).sub (primeSide_tendsto hsuppM hsuppG hpt)
    exact (Complex.continuous_ofReal.tendsto _).comp h
  exact tendsto_nhds_unique (hZ.congr hEF) hW

/-! ## 4. Assembly: `(a) + (b) + (c) → ArchMatch` (PROVED) -/

open GppArchimedeanId

/-- `ψ(1/4 + i r/2)`. -/
noncomputable def dig (r : ℝ) : ℂ := Complex.digamma ((1 / 4 : ℂ) + (r / 2 : ℝ) * Complex.I)

theorem dig_bound : ∃ C : ℝ, 0 ≤ C ∧ ∀ r : ℝ, ‖dig r‖ ≤ C * (1 + |r|) := by
  obtain ⟨C, hC⟩ := GppShiftContour.norm_digamma_strip_le
  refine ⟨max C 0, le_max_right _ _, fun r => ?_⟩
  have h := hC (1 / 2) r ⟨by norm_num, by norm_num⟩
  have e : (((1 / 2 : ℝ) : ℂ) + r * Complex.I) / 2 = (1 / 4 : ℂ) + (r / 2 : ℝ) * Complex.I := by
    push_cast; ring
  rw [e] at h
  exact h.trans (mul_le_mul_of_nonneg_right (le_max_left _ _) (by positivity))

theorem not_neg_nat (r : ℝ) : ∀ m : ℕ, (1 / 4 : ℂ) + (r / 2 : ℝ) * Complex.I ≠ -m := by
  intro m hm
  have h := congrArg Complex.re hm
  simp at h
  have h' : (4 : ℝ) * m = -1 := by linarith
  have : (0 : ℝ) ≤ 4 * m := by positivity
  linarith

theorem dig_continuous : Continuous dig := by
  refine continuous_iff_continuousAt.2 fun r => ?_
  have hz : (1 / 4 : ℂ) + (r / 2 : ℝ) * Complex.I ∈ Complex.integerComplement := by
    rintro ⟨n, hn⟩
    have h := congrArg Complex.re hn
    simp at h
    have h' : (4 * n : ℤ) = 1 := by exact_mod_cast (by linarith : (4 : ℝ) * n = 1)
    omega
  unfold dig
  exact (GppDigammaGaussSeries.analyticAt_digamma hz).continuousAt.comp
    (f := fun r : ℝ => (1 / 4 : ℂ) + ((r / 2 : ℝ) : ℂ) * Complex.I) (by fun_prop)

theorem dig_neg (r : ℝ) : (dig (-r)).re = (dig r).re := by
  unfold dig
  have e : (1 / 4 : ℂ) + ((-r / 2 : ℝ) : ℂ) * Complex.I =
      (starRingEnd ℂ) ((1 / 4 : ℂ) + ((r / 2 : ℝ) : ℂ) * Complex.I) := by
    rw [map_add, map_mul, Complex.conj_ofReal, Complex.conj_I, map_div₀, map_one, map_ofNat]
    push_cast; ring
  rw [e, GppDigammaQuarterSech.digamma_conj (not_neg_nat r), Complex.conj_re]

/-- `ĝ_d(1/2 + ir) = Ψ(r)Ψ(−r) cos(rd)` for real `r ≠ 0` (`integral_gd_exp` + `triFourier`). -/
theorem hline_gd' (d : ℝ) {r : ℝ} (hr : r ≠ 0) :
    hline (gd d) r = efCoeff r * Complex.cos (r * d) := by
  unfold hline GppMeyerSpectralWeil.spectralTransform
  have hz : (1 / 2 : ℂ) + r * Complex.I - 1 / 2 = Complex.I * (r : ℂ) := by ring
  have hr' : (r : ℂ) ≠ 0 := by exact_mod_cast hr
  rw [hz, integral_gd_exp, triFourier _ hr', Complex.cos]
  ring_nf

theorem bigPsi_eq {r : ℝ} (hr : r ≠ 0) : bigPsi r = (Real.sin (r / 2) / (r / 2)) ^ 2 := by
  unfold bigPsi; rw [if_neg hr]

/-- `ĝ_d(1/2 + ir) = Ψ(r)² cos(rd)` for real `r ≠ 0`. -/
theorem hline_gd (d : ℝ) {r : ℝ} (hr : r ≠ 0) :
    hline (gd d) r = ((bigPsi r ^ 2 * Real.cos (r * d) : ℝ) : ℂ) := by
  rw [hline_gd' d hr, bigPsi_eq hr]
  unfold efCoeff GppKreinLaplace.triPsi
  push_cast
  rw [neg_div, Complex.sin_neg, neg_div_neg_eq]
  ring

theorem bigPsi_neg (r : ℝ) : bigPsi (-r) = bigPsi r := by
  by_cases hr : r = 0
  · rw [hr, neg_zero]
  · rw [bigPsi_eq hr, bigPsi_eq (neg_ne_zero.2 hr), neg_div, Real.sin_neg, neg_div_neg_eq]

theorem bigPsi_measurable : Measurable bigPsi := by
  unfold bigPsi
  exact Measurable.ite (measurableSet_singleton (0 : ℝ)) measurable_const (by fun_prop)

theorem norm_hline_gd_le (d : ℝ) {r : ℝ} (hr : r ≠ 0) :
    ‖hline (gd d) r‖ ≤ 144 ^ 2 / (1 + |r|) ^ 4 := by
  rw [hline_gd' d hr, norm_mul]
  have h1 := norm_triPsi_le (γ := (r : ℂ)) (by norm_num)
  have h2 := norm_triPsi_le (γ := -(r : ℂ)) (by norm_num)
  rw [Complex.norm_real, Real.norm_eq_abs] at h1
  rw [norm_neg, Complex.norm_real, Real.norm_eq_abs] at h2
  have hc : ‖Complex.cos (r * d)‖ ≤ 1 := by
    rw [← Complex.ofReal_mul, ← Complex.ofReal_cos, Complex.norm_real, Real.norm_eq_abs]
    exact Real.abs_cos_le_one _
  unfold efCoeff
  rw [norm_mul]
  calc ‖GppKreinLaplace.triPsi r‖ * ‖GppKreinLaplace.triPsi (-r)‖ * ‖Complex.cos (r * d)‖
      ≤ (144 / (1 + |r|) ^ 2) * (144 / (1 + |r|) ^ 2) * 1 := by gcongr
    _ = 144 ^ 2 / (1 + |r|) ^ 4 := by field_simp

theorem ae_ne_zero : ∀ᵐ r : ℝ, r ≠ 0 := compl_mem_ae_iff.2 (measure_singleton 0)

theorem bound_pt (d : ℝ) {w : ℝ → ℂ} {K : ℝ} (_hK : 0 ≤ K) (hw : ∀ r, ‖w r‖ ≤ K * (1 + |r|))
    {r : ℝ} (hr : r ≠ 0) :
    ‖hline (gd d) r‖ * ‖w r‖ ≤ 144 ^ 2 * K * (1 + ‖r‖) ^ (-(3 : ℝ)) := by
  have ha : 0 < 1 + |r| := by positivity
  have hp : (1 + ‖r‖) ^ (-(3 : ℝ)) = 1 / (1 + |r|) ^ 3 := by
    rw [Real.rpow_neg (by positivity), Real.norm_eq_abs, one_div]
    norm_cast
  rw [hp]
  calc ‖hline (gd d) r‖ * ‖w r‖ ≤ (144 ^ 2 / (1 + |r|) ^ 4) * (K * (1 + |r|)) :=
        mul_le_mul (norm_hline_gd_le d hr) (hw r) (norm_nonneg _) (by positivity)
    _ = 144 ^ 2 * K * (1 / (1 + |r|) ^ 3) := by field_simp

theorem hline_gd_aesm (d : ℝ) : AEStronglyMeasurable (hline (gd d)) := by
  have hm : Measurable (fun r : ℝ => ((bigPsi r ^ 2 * Real.cos (r * d) : ℝ) : ℂ)) :=
    Complex.measurable_ofReal.comp ((bigPsi_measurable.pow_const 2).mul (by fun_prop))
  refine hm.aestronglyMeasurable.congr ?_
  filter_upwards [ae_ne_zero] with r hr using (hline_gd d hr).symm

theorem integrable_hline_gd_mul (d : ℝ) {w : ℝ → ℂ} {K : ℝ} (hK : 0 ≤ K)
    (hwc : Continuous w) (hw : ∀ r, ‖w r‖ ≤ K * (1 + |r|)) :
    Integrable (fun r => hline (gd d) r * w r) := by
  refine ((integrable_one_add_norm (E := ℝ) (r := 3) (by norm_num)).const_mul
    (144 ^ 2 * K)).mono' ((hline_gd_aesm d).mul hwc.aestronglyMeasurable) ?_
  filter_upwards [ae_ne_zero] with r hr
  rw [norm_mul]
  exact bound_pt d hK hw hr

/-- Dominated convergence on the critical line, from `(b)`. PROVED (conditional on `(b)`). -/
theorem tendsto_hline_mul (hb : StripFactor) (d : ℝ) {w : ℝ → ℂ} {K : ℝ} (hK : 0 ≤ K)
    (hwc : Continuous w) (hw : ∀ r, ‖w r‖ ≤ K * (1 + |r|)) :
    Tendsto (fun n => ∫ r, hline (moll d n) r * w r) atTop (𝓝 (∫ r, hline (gd d) r * w r)) := by
  obtain ⟨hfac, hbd, hlim⟩ := hb
  have hs : ∀ r : ℝ, 0 ≤ ((1 / 2 : ℂ) + r * Complex.I).re ∧ ((1 / 2 : ℂ) + r * Complex.I).re ≤ 1 :=
    fun r => by simp; norm_num
  have hf : ∀ n r, hline (moll d n) r = hline (gd d) r * mhat n ((1 / 2 : ℂ) + r * Complex.I) :=
    fun n r => hfac d n _ (hs r).1 (hs r).2
  refine tendsto_integral_of_dominated_convergence
    (fun r => Real.exp (1 / 2) * (144 ^ 2 * K * (1 + ‖r‖) ^ (-(3 : ℝ)))) (fun n => ?_)
    ((integrable_one_add_norm (E := ℝ) (r := 3) (by norm_num)).const_mul (144 ^ 2 * K)
      |>.const_mul _) (fun n => ?_) (Eventually.of_forall fun r => ?_)
  · exact ((hline_continuous (moll_isTest d n)).mul hwc).aestronglyMeasurable
  · filter_upwards [ae_ne_zero] with r hr
    rw [hf, norm_mul, norm_mul]
    calc ‖hline (gd d) r‖ * ‖mhat n ((1 / 2 : ℂ) + r * Complex.I)‖ * ‖w r‖
        = ‖mhat n ((1 / 2 : ℂ) + r * Complex.I)‖ * (‖hline (gd d) r‖ * ‖w r‖) := by ring
      _ ≤ Real.exp (1 / 2) * (144 ^ 2 * K * (1 + ‖r‖) ^ (-(3 : ℝ))) :=
        mul_le_mul (hbd n _ (hs r).1 (hs r).2) (bound_pt d hK hw hr) (by positivity)
          (by positivity)
  · have e : ∀ n, hline (moll d n) r * w r =
        (hline (gd d) r * w r) * mhat n ((1 / 2 : ℂ) + r * Complex.I) := fun n => by
      rw [hf]; ring
    simp_rw [e]
    simpa using (hlim _ (hs r).1 (hs r).2).const_mul (hline (gd d) r * w r)

/-- The in-tree `archimedean_identification` extended to `g_d`. PROVED (conditional on (a)–(c)). -/
theorem archId_gd (ha : MollC1) (hb : StripFactor) (hc : RealPlaceDCT) (d : ℝ) :
    ((1 / (2 * Real.pi) : ℝ) : ℂ) * ∫ r : ℝ, hline (gd d) r * dig r =
      ((GppArchimedeanFloor.archimedeanTerm (gd d) + gd d 0 * Real.log Real.pi : ℝ) : ℂ) := by
  obtain ⟨C, hC0, hC⟩ := dig_bound
  have hL := (tendsto_hline_mul hb d hC0 dig_continuous hC).const_mul
    ((1 / (2 * Real.pi) : ℝ) : ℂ)
  have hR : Tendsto (fun n => ((GppArchimedeanFloor.archimedeanTerm (moll d n) +
      moll d n 0 * Real.log Real.pi : ℝ) : ℂ)) atTop
      (𝓝 ((GppArchimedeanFloor.archimedeanTerm (gd d) + gd d 0 * Real.log Real.pi : ℝ) : ℂ)) :=
    (Complex.continuous_ofReal.tendsto _).comp
      ((archLimit ha hc d).add (((ha d).1.tendsto_at 0).mul_const _))
  exact tendsto_nhds_unique (hL.congr fun n => archimedean_identification (moll_isTest d n)) hR

/-- Fourier inversion at `0` for `g_d`, by the same limit. PROVED (conditional on (a), (b)). -/
theorem integral_hline_gd (ha : MollC1) (hb : StripFactor) (d : ℝ) :
    ∫ r : ℝ, hline (gd d) r = ((2 * Real.pi : ℝ) : ℂ) * (gd d 0 : ℂ) := by
  have hL := tendsto_hline_mul hb d (w := fun _ => (1 : ℂ)) zero_le_one continuous_const
    (fun r => by rw [norm_one, one_mul]; linarith [abs_nonneg r])
  simp only [mul_one] at hL
  have hR : Tendsto (fun n => ((2 * Real.pi : ℝ) : ℂ) * (moll d n 0 : ℂ)) atTop
      (𝓝 (((2 * Real.pi : ℝ) : ℂ) * (gd d 0 : ℂ))) :=
    ((Complex.continuous_ofReal.tendsto _).comp ((ha d).1.tendsto_at 0)).const_mul _
  exact tendsto_nhds_unique (hL.congr fun n => integral_hline (moll_isTest d n)) hR

/-- The real frequency integrand `Ψ(r)² cos(rd) (Re ψ(1/4 + ir/2) − log π)`. -/
noncomputable def archF (d r : ℝ) : ℝ :=
  bigPsi r ^ 2 * Real.cos (r * d) * ((dig r).re - Real.log Real.pi)

theorem archF_abs (d r : ℝ) : archF d |r| = archF d r := by
  rcases abs_choice r with h | h
  · rw [h]
  · rw [h]; unfold archF; rw [bigPsi_neg, neg_mul, Real.cos_neg, dig_neg]

theorem archMatch_of_moll (ha : MollC1) (hb : StripFactor) (hc : RealPlaceDCT) : ArchMatch := by
  intro d
  obtain ⟨C, hC0, hC⟩ := dig_bound
  have hI1 := integrable_hline_gd_mul d hC0 dig_continuous hC
  have hI0 : Integrable (hline (gd d)) := by
    have := integrable_hline_gd_mul d (w := fun _ => (1 : ℂ)) zero_le_one continuous_const
      (fun r => by rw [norm_one, one_mul]; linarith [abs_nonneg r])
    simpa only [mul_one] using this
  have hA := archId_gd ha hb hc d
  have h0 := integral_hline_gd ha hb d
  have hpi : (0 : ℝ) < Real.pi := Real.pi_pos
  -- `∫ archF = 2π · archimedeanTerm g_d`
  have hF : ∫ r, archF d r =
      (∫ r, hline (gd d) r * dig r).re - Real.log Real.pi * (∫ r, hline (gd d) r).re := by
    have e1 := integral_re hI1
    have e0 := integral_re hI0
    simp only [RCLike.re_to_complex] at e1 e0
    rw [← e1, ← e0, ← integral_const_mul,
      ← integral_sub (f := fun x => (hline (gd d) x * dig x).re)
        (g := fun a => Real.log Real.pi * (hline (gd d) a).re) hI1.re (hI0.re.const_mul _)]
    refine integral_congr_ae ?_
    filter_upwards [ae_ne_zero] with r hr
    unfold archF
    rw [hline_gd d hr, Complex.re_ofReal_mul, Complex.ofReal_re]
    ring
  have hval : (∫ r, hline (gd d) r * dig r) =
      ((2 * Real.pi : ℝ) : ℂ) *
        ((GppArchimedeanFloor.archimedeanTerm (gd d) + gd d 0 * Real.log Real.pi : ℝ) : ℂ) := by
    rw [← hA]
    have : ((2 * Real.pi : ℝ) : ℂ) * ((1 / (2 * Real.pi) : ℝ) : ℂ) = 1 := by
      rw [← Complex.ofReal_mul]; field_simp; simp
    rw [← mul_assoc, this, one_mul]
  rw [hval, h0, ← Complex.ofReal_mul, ← Complex.ofReal_mul, Complex.ofReal_re,
    Complex.ofReal_re] at hF
  have hfold : ∫ r, archF d r = 2 * ∫ r in Set.Ioi (0 : ℝ), archF d r := by
    rw [← integral_comp_abs]
    exact integral_congr_ae (Eventually.of_forall fun r => (archF_abs d r).symm)
  unfold archTerm
  change _ = (1 / Real.pi) * ∫ t in Set.Ioi (0 : ℝ), archF d t
  have : ∫ t in Set.Ioi (0 : ℝ), archF d t = Real.pi * GppArchimedeanFloor.archimedeanTerm (gd d) := by
    linarith
  rw [this]; field_simp

/-- **`HEFIdent` from the three mollification lemmas (PROVED, conditional).** -/
theorem hefIdent_of_moll (ha : MollC1) (hb : StripFactor) (hc : RealPlaceDCT) : HEFIdent :=
  hefIdent_of_reg_arch (regularity_of_moll ha hb hc) (archMatch_of_moll ha hb hc)

/-- The fine chain with `Regularity` and `ArchMatch` both reduced to (a)–(c). -/
def HPsiMollChain : Prop :=
  GppKreinHPsi.HPsi → MollC1 → StripFactor → RealPlaceDCT → KreinDecomp →
    GppSpectralMatchingSplit.IdentityStep → RiemannHypothesis

theorem hPsiMollChain : HPsiMollChain := fun h ha hb hc hD hI =>
  hPsiRegArchChain h (regularity_of_moll ha hb hc) (archMatch_of_moll ha hb hc) hD hI

end GppMollifierSplit

#print axioms GppMollifierSplit.moll_apply
#print axioms GppMollifierSplit.bump_eq_zero
#print axioms GppMollifierSplit.gd_eq_zero
#print axioms GppMollifierSplit.gd_measurable
#print axioms GppMollifierSplit.abs_gd_le
#print axioms GppMollifierSplit.gd_integrable
#print axioms GppMollifierSplit.moll_eq_zero
#print axioms GppMollifierSplit.moll_isTest
#print axioms GppMollifierSplit.moll_even
#print axioms GppMollifierSplit.abs_moll_le
#print axioms GppMollifierSplit.regSumm
#print axioms GppMollifierSplit.primeSide_tendsto
#print axioms GppMollifierSplit.poleTerm_tendsto
#print axioms GppMollifierSplit.archLimit
#print axioms GppMollifierSplit.regularity_of_moll
#print axioms GppMollifierSplit.dig_bound
#print axioms GppMollifierSplit.not_neg_nat
#print axioms GppMollifierSplit.dig_continuous
#print axioms GppMollifierSplit.dig_neg
#print axioms GppMollifierSplit.hline_gd'
#print axioms GppMollifierSplit.bigPsi_eq
#print axioms GppMollifierSplit.hline_gd
#print axioms GppMollifierSplit.bigPsi_neg
#print axioms GppMollifierSplit.bigPsi_measurable
#print axioms GppMollifierSplit.norm_hline_gd_le
#print axioms GppMollifierSplit.ae_ne_zero
#print axioms GppMollifierSplit.bound_pt
#print axioms GppMollifierSplit.hline_gd_aesm
#print axioms GppMollifierSplit.integrable_hline_gd_mul
#print axioms GppMollifierSplit.tendsto_hline_mul
#print axioms GppMollifierSplit.archId_gd
#print axioms GppMollifierSplit.integral_hline_gd
#print axioms GppMollifierSplit.archF_abs
#print axioms GppMollifierSplit.archMatch_of_moll
#print axioms GppMollifierSplit.hefIdent_of_moll
#print axioms GppMollifierSplit.hPsiMollChain
