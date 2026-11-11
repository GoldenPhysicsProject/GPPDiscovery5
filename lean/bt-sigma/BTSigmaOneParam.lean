import GppVerify.RiemannHypothesis.BTSigmaNonCov

/-!
# BT-Σ: the one-parameter-group test

Scratch module (Heart of the Cosmos / Grok, GPPDiscovery5 seat) for Turbo's question on
GPP-bridge PR #11, comment 6804422975.

**Part 1 (bounded generators, in particular finite-dimensional `Gp`).** `analytic_orbit_kill`:
for independent intertwiners `U_W, V_W` and real characters with `χ₀χ₁ ≠ 0` on `W ≥ 0`, if both
boundary orbits `W ↦ U_W a`, `W ↦ V_W b` are real-analytic on `ℝ`, BT3 is contradictory on any
nonnegative nonzero smooth compactly supported test with a prime `p ∈ S`, `p ≤ N`. The exact step:
BT3 gives `stratumSN(W) = −2χ₀χ₁(W)·Im⟪U_W a, V_W b⟫`; the right factor is real-analytic and
vanishes on `W > max 0 (R + log N)`; the identity theorem
(`AnalyticOnNhd.eqOn_zero_of_preconnected_of_eventuallyEq_zero`) makes it vanish on `ℝ`, so
`stratumSN(log p) = 0`, against `stratumSN_pos`. `oneParam_bounded_kill` /
`oneParam_finiteDim_kill`: `U_W = exp(W•X)`, `V_W = exp(W•Y)` with bounded `X, Y` have analytic
orbits (`analytic_exp_orbit`), so they cannot evade. No analyticity or group law is needed for the
characters.

**Part 2 (infinite-dimensional `Gp`).** `analytic_orbit_kill` holds for any `Gp`, so the kill
extends exactly as far as the boundary orbits are real-analytic (e.g. analytic vectors of the
generators). Without that it does not extend: see the bridge note for the translation-group
evader on `L²(ℝ)` (paper level, not formalised here).

All PROVED, standard axioms. No `sorry`, no `admit`, no custom axiom. S_b existence, (P-COV),
(P-UNIT), (P-LOG), (BT-Σ), ARP and RH remain OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo.
-/

namespace GppBTSigmaOneParam

open GppBTSigmaLocalEuler GppBTSigma GppArithmeticReflection Complex MeasureTheory
open GppBTSigmaCovObstruction GppBTSigmaBoundaryPW GppBTSigmaNonCov
open scoped ContDiff Topology

variable {H Gu Gp : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
  [NormedAddCommGroup Gu] [InnerProductSpace ℂ Gu]
  [NormedAddCommGroup Gp] [InnerProductSpace ℂ Gp]

/-- On the diagonal the block form is `−2i·Im⟪x, y⟫`. PROVED. -/
theorem blockForm_diag (x y : Gp) :
    blockForm x y x y = -2 * I * (((inner ℂ x y).im : ℝ) : ℂ) := by
  unfold blockForm
  have h1 : (inner ℂ y x).re = (inner ℂ x y).re := by
    rw [← inner_conj_symm x y, Complex.conj_re]
  have h2 : (inner ℂ y x).im = -(inner ℂ x y).im := by
    rw [← inner_conj_symm x y, Complex.conj_im, neg_neg]
  apply Complex.ext
  · simp only [sub_re, mul_re, ofReal_re, ofReal_im, I_re, I_im]
    rw [h1]; norm_num
  · simp only [sub_im, mul_im, mul_re, ofReal_re, ofReal_im, I_re, I_im]
    rw [h2]; norm_num; ring

/-- `W ↦ Im⟪x W, y W⟫` is real-analytic when both orbits are. PROVED. -/
theorem analytic_im_inner {x y : ℝ → Gp} (hx : AnalyticOnNhd ℝ x Set.univ)
    (hy : AnalyticOnNhd ℝ y Set.univ) :
    AnalyticOnNhd ℝ (fun W => (inner ℂ (x W) (y W)).im) Set.univ := by
  let _ : InnerProductSpace ℝ Gp := InnerProductSpace.complexToReal
  have e : (fun W => (inner ℂ (x W) (y W)).im) =
      fun W => (innerSL ℝ (x W)) ((-I : ℂ) • y W) := by
    funext W
    rw [innerSL_apply_apply, real_inner_eq_re_inner ℂ, inner_smul_right]
    simp
  rw [e]
  intro W _
  have hp : AnalyticAt ℝ (fun W => (x W, (-I : ℂ) • y W)) W :=
    AnalyticAt.prod (hx W trivial) (AnalyticAt.const_smul (hy W trivial))
  exact (ContinuousLinearMap.analyticAt_bilinear (innerSL ℝ : Gp →L[ℝ] Gp →L[ℝ] ℝ) _).comp hp

/-- **The one-parameter-group kill (analytic orbits).** Independent intertwiners `U_W`, `V_W` and
any real characters with `χ₀χ₁ ≠ 0` on `W ≥ 0`. If both boundary orbits `W ↦ U_W a`,
`W ↦ V_W b` are real-analytic on `ℝ`, then BT3 plus a nonnegative nonzero test with a prime
`p ∈ S`, `p ≤ N` is contradictory. The step: `Im⟪U_W a, V_W b⟫` is real-analytic and vanishes on
`W > max 0 (R + log N)`, so by the identity theorem it vanishes on `ℝ`, which contradicts
`stratumSN_pos` at `W = log p`. PROVED. -/
theorem analytic_orbit_kill (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S) (hpN : p ≤ N)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hf0 : ∀ x, 0 ≤ f x) {x₀ : ℝ} (hx : f x₀ ≠ 0)
    (U V : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ)
    (hcov : OrbitIndepCovariantPrime D (autocorr f) U V χ₀ χ₁)
    (hχ : ∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0)
    (hU : AnalyticOnNhd ℝ (fun W => U W (D.Γ₀p (D.testVec (autocorr f) 0))) Set.univ)
    (hV : AnalyticOnNhd ℝ (fun W => V W (D.Γ₁p (D.testVec (autocorr f) 0))) Set.univ) :
    False := by
  set a := D.Γ₀p (D.testVec (autocorr f) 0)
  set b := D.Γ₁p (D.testVec (autocorr f) 0)
  set g : ℝ → ℝ := fun W => (inner ℂ (U W a) (V W b)).im with hgdef
  have hg : AnalyticOnNhd ℝ g Set.univ := analytic_im_inner hU hV
  have key : ∀ W : ℝ, 0 ≤ W →
      stratumSN S N (autocorr f) W = -2 * (χ₀ W * χ₁ W) * g W := by
    intro W hW
    have e := hP f hf hc W hW
    rw [primeBlock_indep D _ U V χ₀ χ₁ hcov W hW, blockForm_diag] at e
    have e2 : ((stratumSN S N (autocorr f) W : ℝ) : ℂ) =
        ((-2 * (χ₀ W * χ₁ W) * g W : ℝ) : ℂ) := by
      apply mul_left_cancel₀ Complex.I_ne_zero
      rw [← e]; simp only [hgdef]; push_cast; ring
    exact_mod_cast e2
  obtain ⟨R, hR⟩ := autocorr_eventually_zero hc
  set T : ℝ := max 0 (R + Real.log N)
  have htail : ∀ W : ℝ, T < W → g W = 0 := by
    intro W hW
    have hW0 : 0 ≤ W := le_trans (le_max_left _ _) hW.le
    have hWR : R + Real.log N < W := lt_of_le_of_lt (le_max_right _ _) hW
    have e := key W hW0
    rw [stratumSN_eq_zero_of_large S N _ R hR W hWR] at e
    rcases mul_eq_zero.mp e.symm with h1 | h1
    · exact absurd (by linarith : χ₀ W * χ₁ W = 0) (hχ W hW0)
    · exact h1
  have hev : g =ᶠ[𝓝 (T + 1)] 0 := by
    filter_upwards [Ioi_mem_nhds (by linarith : T < T + 1)] with W hW
    exact htail W hW
  have hz := hg.eqOn_zero_of_preconnected_of_eventuallyEq_zero isPreconnected_univ
    (Set.mem_univ _) hev
  have hlog := key (Real.log p) (Real.log_natCast_nonneg p)
  rw [show g (Real.log p) = 0 from hz (Set.mem_univ _), mul_zero] at hlog
  have hpos := stratumSN_pos S N hp hpS hpN hf.continuous hc hf0 hx
  linarith

/-- The orbit `W ↦ exp(W•X) v` of a bounded generator is real-analytic on `ℝ`. PROVED. -/
theorem analytic_exp_orbit [CompleteSpace Gp] (X : Gp →L[ℂ] Gp) (v : Gp) :
    AnalyticOnNhd ℝ (fun W : ℝ => NormedSpace.exp ((W : ℂ) • X) v) Set.univ := by
  intro W _
  have h1 : AnalyticAt ℂ (fun z : ℂ => NormedSpace.exp (z • X) v) (W : ℂ) := by
    have hs : AnalyticAt ℂ (fun z : ℂ => z • X) (W : ℂ) :=
      ((ContinuousLinearMap.id ℂ ℂ).smulRight X).analyticAt _
    have he : AnalyticAt ℂ (fun z : ℂ => NormedSpace.exp (z • X)) (W : ℂ) :=
      AnalyticAt.comp (f := fun z : ℂ => z • X) (NormedSpace.exp_analytic (𝕂 := ℂ) _) hs
    exact AnalyticAt.comp (f := fun z : ℂ => NormedSpace.exp (z • X))
      ((ContinuousLinearMap.apply ℂ Gp v).analyticAt _) he
  have h2 : AnalyticAt ℝ (fun z : ℂ => NormedSpace.exp (z • X) v) (W : ℂ) :=
    h1.restrictScalars
  exact h2.comp (Complex.ofRealCLM.analyticAt W)

/-- **Part 1: bounded generators cannot evade.** If `U_W = exp(W•X)` and `V_W = exp(W•Y)` with
bounded generators `X, Y` (every one-parameter group on a finite-dimensional `Gp` has this form),
then under BT3 no real characters with `χ₀χ₁ ≠ 0` make `Γ₀p, Γ₁p` independently covariant along a
nonnegative nonzero test orbit. The relative pairing is entire in `W` (a trigonometric polynomial
when `Gp` is finite-dimensional and `X, Y` skew-adjoint). PROVED. -/
theorem oneParam_bounded_kill [CompleteSpace Gp] (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ)
    (N : ℕ) (hP : BTSigmaTraceSN D S N) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S) (hpN : p ≤ N)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hf0 : ∀ x, 0 ≤ f x) {x₀ : ℝ} (hx : f x₀ ≠ 0)
    (U V : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ)
    (hcov : OrbitIndepCovariantPrime D (autocorr f) U V χ₀ χ₁)
    (hχ : ∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0) (X Y : Gp →L[ℂ] Gp)
    (hUe : ∀ (W : ℝ) (v : Gp), U W v = NormedSpace.exp ((W : ℂ) • X) v)
    (hVe : ∀ (W : ℝ) (v : Gp), V W v = NormedSpace.exp ((W : ℂ) • Y) v) : False := by
  refine analytic_orbit_kill D S N hP hp hpS hpN f hf hc hf0 hx U V χ₀ χ₁ hcov hχ ?_ ?_
  · simpa only [hUe] using analytic_exp_orbit X (D.Γ₀p (D.testVec (autocorr f) 0))
  · simpa only [hVe] using analytic_exp_orbit Y (D.Γ₁p (D.testVec (autocorr f) 0))

/-- Finite-dimensional `Gp` is complete, so `oneParam_bounded_kill` applies without a separate
completeness hypothesis. PROVED. -/
theorem oneParam_finiteDim_kill [FiniteDimensional ℂ Gp] (D : BTSigmaDatum H Gu Gp)
    (S : Finset ℕ) (N : ℕ) (hP : BTSigmaTraceSN D S N) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S)
    (hpN : p ≤ N) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hf0 : ∀ x, 0 ≤ f x) {x₀ : ℝ} (hx : f x₀ ≠ 0)
    (U V : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ)
    (hcov : OrbitIndepCovariantPrime D (autocorr f) U V χ₀ χ₁)
    (hχ : ∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0) (X Y : Gp →L[ℂ] Gp)
    (hUe : ∀ (W : ℝ) (v : Gp), U W v = NormedSpace.exp ((W : ℂ) • X) v)
    (hVe : ∀ (W : ℝ) (v : Gp), V W v = NormedSpace.exp ((W : ℂ) • Y) v) : False :=
  haveI : CompleteSpace Gp := FiniteDimensional.complete ℂ Gp
  oneParam_bounded_kill D S N hP hp hpS hpN f hf hc hf0 hx U V χ₀ χ₁ hcov hχ X Y hUe hVe

end GppBTSigmaOneParam

#print axioms GppBTSigmaOneParam.blockForm_diag
#print axioms GppBTSigmaOneParam.analytic_im_inner
#print axioms GppBTSigmaOneParam.analytic_orbit_kill
#print axioms GppBTSigmaOneParam.analytic_exp_orbit
#print axioms GppBTSigmaOneParam.oneParam_bounded_kill
#print axioms GppBTSigmaOneParam.oneParam_finiteDim_kill
