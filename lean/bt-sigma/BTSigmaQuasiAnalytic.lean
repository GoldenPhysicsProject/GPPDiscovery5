import GppVerify.RiemannHypothesis.BTSigmaOneParam

/-!
# BT-Σ: the quasi-analytic (Denjoy–Carleman) BT3 test

Scratch module (Heart of the Cosmos / Grok, GPPDiscovery5 seat) for Turbo's item 3 on
GPP-bridge PR #11 (comments 6806968693, 6807522996): state the exact Denjoy–Carleman class and
test BT3 against it.

Setting (as in `BTSigmaOneParam`): independent intertwiners `U_W, V_W` on `Gp`, real characters
with `χ₀χ₁ ≠ 0` on `W ≥ 0`, boundary vectors `a = Γ₀p ψ₀`, `b = Γ₁p ψ₀` for the test
`ψ₀ = testVec (autocorr f) 0`, and the relative pairing `g(W) = Im⟪U_W a, V_W b⟫`.

**The class.** For a weight `M : ℕ → ℝ` and a set `I`, `dcClass M I` is the Denjoy–Carleman class
`C{M}(I)`: smooth `g` with `|g⁽ⁿ⁾(x)| ≤ C ρⁿ Mₙ` for all `n` and `x ∈ I`. The weight is
*quasi-analytic* (`DCQuasiAnalyticWeight`) when `M₀ = 1`, `Mₙ > 0`, `M` is log-convex and
`∑ Mₙ/Mₙ₊₁ = ∞`. The Denjoy–Carleman theorem says that then `C{M}(I)` has the uniqueness property
on an interval (`QuasiAnalyticOn`). That theorem is cited, not formalised: it is the Prop
`DenjoyCarlemanUniqueness`, used only as a hypothesis.

**Results (all PROVED, standard axioms).**
* `bt3_pairing_tail`: BT3 forces `g = 0` on `W > T` and `g(log p) ≠ 0`. This is the exact step.
* `quasiAnalytic_orbit_kill`: for ANY class with the uniqueness property on `[0, ∞)`, BT3 is
  contradictory if `g` lies in the class.
* `bt3_forces_not_quasiAnalytic`: equivalently, BT3 forces `g` OUT of every quasi-analytic class.
  So "BT3 forces the pairing into a Denjoy–Carleman class" is KILLED.
* `dc_orbit_kill`: the DC corollary, conditional on the cited DC uniqueness for the weight.
* `analytic_quasiAnalyticOn`: real-analytic functions on `ℝ` have the uniqueness property on
  `[0, ∞)` (identity theorem), so the abstract kill recovers `analytic_orbit_kill`.
* `factorial_quasiAnalyticWeight`: `Mₙ = n!` (the analytic class) is a quasi-analytic weight.
* `factorial_sq_not_quasiAnalytic`: `Mₙ = (n!)²` (Gevrey-2) fails the divergence condition, so
  it is non-quasi-analytic, consistent with compactly supported profiles living there.

No `sorry`, no `admit`, no custom axiom. S_b existence, (P-COV), (P-UNIT), (P-LOG), (BT-Σ), ARP
and RH remain OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo.
-/

namespace GppBTSigmaQuasiAnalytic

open GppBTSigmaLocalEuler GppBTSigma GppArithmeticReflection Complex MeasureTheory
open GppBTSigmaCovObstruction GppBTSigmaBoundaryPW GppBTSigmaNonCov GppBTSigmaOneParam
open scoped ContDiff Topology Nat

variable {H Gu Gp : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
  [NormedAddCommGroup Gu] [InnerProductSpace ℂ Gu]
  [NormedAddCommGroup Gp] [InnerProductSpace ℂ Gp]

/-- Uniqueness (quasi-analytic) property of a class `𝒞` on `I`: a member that vanishes near an
interior point of `I` vanishes on all of `I`. -/
def QuasiAnalyticOn (𝒞 : Set (ℝ → ℝ)) (I : Set ℝ) : Prop :=
  ∀ g ∈ 𝒞, ∀ x₀ ∈ interior I, g =ᶠ[𝓝 x₀] 0 → Set.EqOn g 0 I

/-- The Denjoy–Carleman class `C{M}(I)`. -/
def dcClass (M : ℕ → ℝ) (I : Set ℝ) : Set (ℝ → ℝ) :=
  {g | ContDiff ℝ ∞ g ∧ ∃ C ρ : ℝ, ∀ n : ℕ, ∀ x ∈ I, ‖iteratedDeriv n g x‖ ≤ C * ρ ^ n * M n}

/-- Quasi-analytic Denjoy–Carleman weight: normalised, positive, log-convex, `∑ Mₙ/Mₙ₊₁ = ∞`. -/
def DCQuasiAnalyticWeight (M : ℕ → ℝ) : Prop :=
  M 0 = 1 ∧ (∀ n, 0 < M n) ∧ (∀ n, M (n + 1) ^ 2 ≤ M n * M (n + 2)) ∧
    ¬ Summable (fun n => M n / M (n + 1))

/-- The Denjoy–Carleman theorem (uniqueness half) for the weight `M` on `I`. Cited, NOT
formalised here; used only as a hypothesis. -/
def DenjoyCarlemanUniqueness (M : ℕ → ℝ) (I : Set ℝ) : Prop :=
  DCQuasiAnalyticWeight M → QuasiAnalyticOn (dcClass M I) I

/-- **The exact step.** BT3 with independent intertwiners forces the relative pairing to vanish
on the tail `W > T` while `g(log p) ≠ 0`. PROVED. -/
theorem bt3_pairing_tail (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S) (hpN : p ≤ N)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hf0 : ∀ x, 0 ≤ f x) {x₀ : ℝ} (hx : f x₀ ≠ 0)
    (U V : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ)
    (hcov : OrbitIndepCovariantPrime D (autocorr f) U V χ₀ χ₁)
    (hχ : ∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0) :
    ∃ T : ℝ, 0 ≤ T ∧
      (∀ W : ℝ, T < W → (inner ℂ (U W (D.Γ₀p (D.testVec (autocorr f) 0)))
        (V W (D.Γ₁p (D.testVec (autocorr f) 0)))).im = 0) ∧
      (inner ℂ (U (Real.log p) (D.Γ₀p (D.testVec (autocorr f) 0)))
        (V (Real.log p) (D.Γ₁p (D.testVec (autocorr f) 0)))).im ≠ 0 := by
  set a := D.Γ₀p (D.testVec (autocorr f) 0)
  set b := D.Γ₁p (D.testVec (autocorr f) 0)
  set g : ℝ → ℝ := fun W => (inner ℂ (U W a) (V W b)).im with hgdef
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
  refine ⟨max 0 (R + Real.log N), le_max_left _ _, ?_, ?_⟩
  · intro W hW
    have hW0 : 0 ≤ W := le_trans (le_max_left _ _) hW.le
    have hWR : R + Real.log N < W := lt_of_le_of_lt (le_max_right _ _) hW
    have e := key W hW0
    rw [stratumSN_eq_zero_of_large S N _ R hR W hWR] at e
    rcases mul_eq_zero.mp e.symm with h1 | h1
    · exact absurd (by linarith : χ₀ W * χ₁ W = 0) (hχ W hW0)
    · exact h1
  · intro h0
    have hlog := key (Real.log p) (Real.log_natCast_nonneg p)
    have h0' : g (Real.log p) = 0 := h0
    rw [h0', mul_zero] at hlog
    have hpos := stratumSN_pos S N hp hpS hpN hf.continuous hc hf0 hx
    linarith

/-- **Quasi-analytic orbit kill.** For any class `𝒞` with the uniqueness property on `[0, ∞)`,
BT3 is contradictory when the relative pairing lies in `𝒞`. PROVED. -/
theorem quasiAnalytic_orbit_kill (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S) (hpN : p ≤ N)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hf0 : ∀ x, 0 ≤ f x) {x₀ : ℝ} (hx : f x₀ ≠ 0)
    (U V : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ)
    (hcov : OrbitIndepCovariantPrime D (autocorr f) U V χ₀ χ₁)
    (hχ : ∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0)
    (𝒞 : Set (ℝ → ℝ)) (h𝒞 : QuasiAnalyticOn 𝒞 (Set.Ici 0))
    (hg : (fun W => (inner ℂ (U W (D.Γ₀p (D.testVec (autocorr f) 0)))
        (V W (D.Γ₁p (D.testVec (autocorr f) 0)))).im) ∈ 𝒞) :
    False := by
  obtain ⟨T, hT0, htail, hlog⟩ :=
    bt3_pairing_tail D S N hP hp hpS hpN f hf hc hf0 hx U V χ₀ χ₁ hcov hχ
  have hint : T + 1 ∈ interior (Set.Ici (0 : ℝ)) := by
    rw [interior_Ici]; exact Set.mem_Ioi.mpr (by linarith)
  have hev : (fun W => (inner ℂ (U W (D.Γ₀p (D.testVec (autocorr f) 0)))
      (V W (D.Γ₁p (D.testVec (autocorr f) 0)))).im) =ᶠ[𝓝 (T + 1)] 0 := by
    filter_upwards [Ioi_mem_nhds (by linarith : T < T + 1)] with W hW
    exact htail W hW
  have hz := h𝒞 _ hg (T + 1) hint hev
  exact hlog (hz (Set.mem_Ici.mpr (Real.log_natCast_nonneg p)))

/-- **BT3 forces the pairing OUT of every quasi-analytic class on `[0, ∞)`.** So the forcing
"BT3 ⇒ pairing ∈ C{M}" is KILLED for every quasi-analytic `M`. PROVED. -/
theorem bt3_forces_not_quasiAnalytic (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S) (hpN : p ≤ N)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hf0 : ∀ x, 0 ≤ f x) {x₀ : ℝ} (hx : f x₀ ≠ 0)
    (U V : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ)
    (hcov : OrbitIndepCovariantPrime D (autocorr f) U V χ₀ χ₁)
    (hχ : ∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0)
    (𝒞 : Set (ℝ → ℝ)) (h𝒞 : QuasiAnalyticOn 𝒞 (Set.Ici 0)) :
    (fun W => (inner ℂ (U W (D.Γ₀p (D.testVec (autocorr f) 0)))
        (V W (D.Γ₁p (D.testVec (autocorr f) 0)))).im) ∉ 𝒞 :=
  fun hg => quasiAnalytic_orbit_kill D S N hP hp hpS hpN f hf hc hf0 hx U V χ₀ χ₁ hcov hχ 𝒞 h𝒞 hg

/-- **Denjoy–Carleman corollary**, conditional on the cited DC uniqueness theorem for `M` on
`[0, ∞)`: a quasi-analytic DC pairing kills BT3. PROVED (given the hypothesis). -/
theorem dc_orbit_kill (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S) (hpN : p ≤ N)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hf0 : ∀ x, 0 ≤ f x) {x₀ : ℝ} (hx : f x₀ ≠ 0)
    (U V : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ)
    (hcov : OrbitIndepCovariantPrime D (autocorr f) U V χ₀ χ₁)
    (hχ : ∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0)
    (M : ℕ → ℝ) (hM : DCQuasiAnalyticWeight M)
    (hDC : DenjoyCarlemanUniqueness M (Set.Ici 0))
    (hg : (fun W => (inner ℂ (U W (D.Γ₀p (D.testVec (autocorr f) 0)))
        (V W (D.Γ₁p (D.testVec (autocorr f) 0)))).im) ∈ dcClass M (Set.Ici 0)) :
    False :=
  quasiAnalytic_orbit_kill D S N hP hp hpS hpN f hf hc hf0 hx U V χ₀ χ₁ hcov hχ _ (hDC hM) hg

/-- Real-analytic functions on `ℝ` have the uniqueness property on `[0, ∞)`. PROVED. -/
theorem analytic_quasiAnalyticOn :
    QuasiAnalyticOn {g : ℝ → ℝ | AnalyticOnNhd ℝ g Set.univ} (Set.Ici 0) := by
  intro g hg x₀ _ hev x _
  exact (AnalyticOnNhd.eqOn_zero_of_preconnected_of_eventuallyEq_zero (f := g) hg
    isPreconnected_univ (Set.mem_univ x₀) hev) (Set.mem_univ x)

/-- The analytic weight `Mₙ = n!` is a quasi-analytic DC weight. PROVED. -/
theorem factorial_quasiAnalyticWeight :
    DCQuasiAnalyticWeight (fun n => (n ! : ℝ)) := by
  refine ⟨by simp, fun n => ?_, fun n => ?_, ?_⟩
  · show (0 : ℝ) < (n ! : ℝ)
    exact_mod_cast Nat.factorial_pos n
  · show ((n + 1)! : ℝ) ^ 2 ≤ (n ! : ℝ) * ((n + 2)! : ℝ)
    have h1 : ((n + 2)! : ℝ) = (n + 2) * (n + 1) * n ! := by
      push_cast [Nat.factorial_succ]; ring
    have h2 : ((n + 1)! : ℝ) = (n + 1) * n ! := by
      push_cast [Nat.factorial_succ]; ring
    rw [h1, h2]
    have hn : (0 : ℝ) ≤ (n ! : ℝ) := by positivity
    nlinarith [mul_nonneg (mul_nonneg hn hn) (by linarith : (0:ℝ) ≤ (n:ℝ) + 1)]
  · have e : (fun n : ℕ => (n ! : ℝ) / ((n + 1)! : ℝ)) = fun n : ℕ => 1 / ((n : ℝ) + 1) := by
      funext n
      have hpos : (0 : ℝ) < (n ! : ℝ) := by exact_mod_cast Nat.factorial_pos n
      rw [Nat.factorial_succ]; push_cast; field_simp
    show ¬ Summable (fun n : ℕ => (n ! : ℝ) / ((n + 1)! : ℝ))
    rw [e]
    intro hs
    apply Real.not_summable_one_div_natCast
    exact (summable_nat_add_iff 1).mp (by simpa [Nat.cast_add, Nat.cast_one] using hs)

/-- The Gevrey-2 weight `Mₙ = (n!)²` is NOT quasi-analytic: `∑ Mₙ/Mₙ₊₁ < ∞`. PROVED. -/
theorem factorial_sq_not_quasiAnalytic :
    ¬ DCQuasiAnalyticWeight (fun n => ((n ! : ℝ)) ^ 2) := by
  rintro ⟨-, -, -, hns⟩
  apply hns
  have e : (fun n : ℕ => (n ! : ℝ) ^ 2 / ((n + 1)! : ℝ) ^ 2) =
      fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 2 := by
    funext n
    have hpos : (0 : ℝ) < (n ! : ℝ) := by exact_mod_cast Nat.factorial_pos n
    rw [Nat.factorial_succ]; push_cast; field_simp
  rw [e]
  have h := (summable_nat_add_iff 1).mpr (Real.summable_one_div_nat_pow.mpr one_lt_two)
  simpa [Nat.cast_add, Nat.cast_one] using h

#print axioms bt3_pairing_tail
#print axioms quasiAnalytic_orbit_kill
#print axioms bt3_forces_not_quasiAnalytic
#print axioms dc_orbit_kill
#print axioms analytic_quasiAnalyticOn
#print axioms factorial_quasiAnalyticWeight
#print axioms factorial_sq_not_quasiAnalytic

end GppBTSigmaQuasiAnalytic
