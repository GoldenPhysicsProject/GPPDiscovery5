import GppVerify.RiemannHypothesis.BTSigmaBoundaryPW
import GppVerify.RiemannHypothesis.BTSigmaCovHalfType

/-!
# BT-Σ: the non-covariant route — two isometries on the prime maps, and the weakened BT1

Scratch module (Heart of the Cosmos / Grok, GPPDiscovery5 seat) for Turbo's instruction on
GPP-bridge PR #11, comment 6804150478 (4).

**Scope.** Fixed test orbit `ψ_W = testVec h W`; no zeros used as input; no existence claim for
`S_b` or for a genuine non-circular datum beyond what is proved below.

## (a) Two isometries `V_W` on `Γ₀p`, `U_W` on `Γ₁p`. Verdict: the obstruction FAILS.
* `primeBlock_twoIsometry` (PROVED): `primeBlock(W) = χ₀χ₁(W)·(⟪U_W b, V_W a⟫ − ⟪V_W a, U_W b⟫)`,
  `a = Γ₀p ψ_0`, `b = Γ₁p ψ_0`.
* `twoIsometry_relPhase_forces_stratum_zero` (PROVED): the BoundaryPW no-go survives under
  `RelPhaseTrivial`: `⟪U_W b, V_W a⟫ = ⟪b, a⟫` on `W ≥ 0`. `relPhase_of_same` (PROVED): this is
  automatic when `V = U` (`LinearIsometry.inner_map_map`). That is the step that breaks: with
  `V ≠ U` the pairing sees the relative operator `V_W* U_W`, not `⟪b, a⟫`.
* `witness_block` (PROVED): on `ℂ²`, `Γ₀p = e₀`, `Γ₁p = (1+|s|)·diag(λ(s),1)e₀` with
  `λ(s) = √(1−t²) + i t`, `t = −s/(2(1+|s|))`, gives `blockForm = I·s` for **any** real `s`.
  `diagIso_nonscalar` (PROVED): `diag(λ,1)` is non-scalar for `λ ≠ 1`.
* `retest_*`, `retest_separates` (PROVED, insertion level): for **any** triplet `D` (its own
  `T*`, `Γ₀, Γ₁`, Green identity, surjectivity) with `Gp = ℂ²` and a unit vector in `Gu`,
  re-choosing only the test vectors by surjectivity gives a datum with finite-cut BT3 and
  (P-UNIT) for every test; two-isometry covariance (`V_W = id`, `U_W = diag(λ_W λ̄_0, 1)`,
  `χ₀ = 1`, `χ₁ > 0`) holds on every orbit; same-isometry covariance fails on a smooth test
  (BoundaryPW); and on that test `V_W* U_W` is non-scalar at some `W ≥ 0`.
  **Caveat (circularity guard):** this is a test-vector insertion, not a non-circular datum;
  (P-COV) is not claimed for the re-test, and the rotation `U_W` depends on the test `h`.

## (b) Weakened BT1. Verdict: HOLDS (nothing banked breaks).
* `BT1Weak` (definition): (P-COV) orbit law + same-isometry intertwining on the unit block +
  two-isometry intertwining (nonvanishing real characters) on the prime block; no `T*`
  covariance.
* `bt1Weak_unitBlock` (PROVED): the unit clause scales the unit block by `α₀α₁`, compatible with
  (P-UNIT) (`punit_unitBlock_covariant`).
* `bt1Weak_growth`, `bt1Weak_type_lt_half`, `bt1Weak_halfType` (PROVED): CovGrowth and
  CovHalfType go through verbatim; they never used the prime intertwining. HalfType's refutation
  of exact half-homogeneity needs no BT1 at all.
* `tstarCov_greenInvariant`, `bt1Weak_obstruction` (PROVED): if `T*`-covariance is kept, the P-COV
  obstruction still fires (`∫ h e^{-v/2} = 0`). So weakening the prime clause alone does not
  evade the obstruction; `T*` must be non-covariant too, i.e. carry the perturbation `P_W` of
  CovPerturbation.
* Consistency: `retest_separates` shows BT3 + (P-UNIT) + the BT1Weak prime clause are jointly
  satisfiable on any `ℂ²` triplet (insertion level); together with the (P-COV) orbit law it is
  OPEN.

S_b existence, a non-circular non-covariant datum, (P-COV), (P-UNIT), (P-LOG), (P-θ), (P-LIM),
(BT-Σ), ARP and RH remain **OPEN**.

No `sorry`, no `admit`, no custom axiom.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo.
-/

namespace GppBTSigmaTwoIsometry

open GppBTSigmaLocalEuler GppBTSigma GppArithmeticReflection Complex MeasureTheory
open GppBTSigmaCovGrowth GppBTSigmaCovObstruction GppBTSigmaCovPerturbation GppBTSigmaCovHalfType
open GppBTSigmaBoundaryPW
open scoped ContDiff InnerProductSpace

section Datum

variable {H Gu Gp : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
  [NormedAddCommGroup Gu] [InnerProductSpace ℂ Gu]
  [NormedAddCommGroup Gp] [InnerProductSpace ℂ Gp]

/-- Two-isometry prime boundary along the orbit: `Γ₀p ψ_W = χ₀(W) V_W Γ₀p ψ_0`,
`Γ₁p ψ_W = χ₁(W) U_W Γ₁p ψ_0`. -/
def TwoIsometryPrime (D : BTSigmaDatum H Gu Gp) (h : ℝ → ℝ) (V U : ℝ → Gp →ₗᵢ[ℂ] Gp)
    (χ₀ χ₁ : ℝ → ℝ) : Prop :=
  ∀ W : ℝ, 0 ≤ W →
    D.Γ₀p (D.testVec h W) = (χ₀ W : ℂ) • V W (D.Γ₀p (D.testVec h 0)) ∧
    D.Γ₁p (D.testVec h W) = (χ₁ W : ℂ) • U W (D.Γ₁p (D.testVec h 0))

/-- The BoundaryPW covariance is the diagonal case `V = U`. PROVED. -/
theorem covariant_iff_twoIsometry_diag (D : BTSigmaDatum H Gu Gp) (h : ℝ → ℝ)
    (U : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ) :
    OrbitCovariantPrime D h U χ₀ χ₁ ↔ TwoIsometryPrime D h U U χ₀ χ₁ := Iff.rfl

/-- The prime block under two isometries. PROVED. -/
theorem primeBlock_twoIsometry (D : BTSigmaDatum H Gu Gp) (h : ℝ → ℝ)
    (V U : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ) (hcov : TwoIsometryPrime D h V U χ₀ χ₁)
    (W : ℝ) (hW : 0 ≤ W) :
    primeBlock D h W = ((χ₀ W * χ₁ W : ℝ) : ℂ) *
      (⟪U W (D.Γ₁p (D.testVec h 0)), V W (D.Γ₀p (D.testVec h 0))⟫_ℂ -
        ⟪V W (D.Γ₀p (D.testVec h 0)), U W (D.Γ₁p (D.testVec h 0))⟫_ℂ) := by
  obtain ⟨e0, e1⟩ := hcov W hW
  unfold primeBlock blockForm
  rw [e0, e1]
  simp only [inner_smul_left, inner_smul_right, Complex.conj_ofReal]
  push_cast; ring

/-- Any orbit on which the prime block scales by a nonvanishing real character has zero stratum
(finite-cut BT3). PROVED. -/
theorem scaling_forces_stratum_zero (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (R : ℝ) (hR : ∀ v, R < v → autocorr f v = 0) (κ : ℝ → ℝ)
    (hscale : ∀ W : ℝ, 0 ≤ W →
      primeBlock D (autocorr f) W = (κ W : ℂ) * primeBlock D (autocorr f) 0)
    (hκ : ∀ W : ℝ, 0 ≤ W → κ W ≠ 0) :
    ∀ W : ℝ, 0 ≤ W → stratumSN S N (autocorr f) W = 0 := by
  have key : ∀ W : ℝ, 0 ≤ W →
      stratumSN S N (autocorr f) W = κ W * stratumSN S N (autocorr f) 0 := by
    intro W hW
    have e := hscale W hW
    rw [hP f hf hc W hW, hP f hf hc 0 le_rfl] at e
    have e2 : ((stratumSN S N (autocorr f) W : ℝ) : ℂ) =
        ((κ W * stratumSN S N (autocorr f) 0 : ℝ) : ℂ) := by
      apply mul_left_cancel₀ Complex.I_ne_zero
      rw [e]; push_cast; ring
    exact_mod_cast e2
  set W₀ : ℝ := max 0 (R + Real.log N + 1)
  have hW₀ : 0 ≤ W₀ := le_max_left _ _
  have hW₀R : R + Real.log N < W₀ := by
    have := le_max_right 0 (R + Real.log N + 1); linarith
  have h0 : stratumSN S N (autocorr f) 0 = 0 := by
    have e := key W₀ hW₀
    rw [stratumSN_eq_zero_of_large S N _ R hR W₀ hW₀R] at e
    rcases mul_eq_zero.mp e.symm with h1 | h1
    · exact absurd h1 (hκ W₀ hW₀)
    · exact h1
  intro W hW
  rw [key W hW, h0, mul_zero]

/-- The relative-phase condition: `V_W* U_W` acts trivially on the pairing of `(a, b)`. -/
def RelPhaseTrivial (D : BTSigmaDatum H Gu Gp) (h : ℝ → ℝ) (V U : ℝ → Gp →ₗᵢ[ℂ] Gp) : Prop :=
  ∀ W : ℝ, 0 ≤ W →
    ⟪U W (D.Γ₁p (D.testVec h 0)), V W (D.Γ₀p (D.testVec h 0))⟫_ℂ =
      ⟪D.Γ₁p (D.testVec h 0), D.Γ₀p (D.testVec h 0)⟫_ℂ

/-- **The obstruction survives when the relative phase is trivial.** PROVED. -/
theorem twoIsometry_relPhase_forces_stratum_zero (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ)
    (N : ℕ) (hP : BTSigmaTraceSN D S N) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f)
    (hc : HasCompactSupport f) (R : ℝ) (hR : ∀ v, R < v → autocorr f v = 0)
    (V U : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ)
    (hcov : TwoIsometryPrime D (autocorr f) V U χ₀ χ₁)
    (hrel : RelPhaseTrivial D (autocorr f) V U)
    (hχ : ∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0) :
    ∀ W : ℝ, 0 ≤ W → stratumSN S N (autocorr f) W = 0 := by
  refine scaling_forces_stratum_zero D S N hP f hf hc R hR (fun W => χ₀ W * χ₁ W) ?_ hχ
  intro W hW
  rw [primeBlock_twoIsometry D _ V U χ₀ χ₁ hcov W hW, hrel W hW,
    ← inner_conj_symm (V W _) (U W _), hrel W hW, inner_conj_symm]
  unfold primeBlock blockForm
  rfl

/-- Same isometry on both maps: the relative phase is automatically trivial. PROVED. -/
theorem relPhase_of_same (D : BTSigmaDatum H Gu Gp) (h : ℝ → ℝ) (U : ℝ → Gp →ₗᵢ[ℂ] Gp) :
    RelPhaseTrivial D h U U := fun _ _ => LinearIsometry.inner_map_map _ _ _

/-! ## (b) The weakened BT1 clause -/

/-- **Weakened BT1 (non-covariant prime block).** (P-COV) orbit law; the unit block intertwines
with one isometry `Uu_W` and real characters `α₀, α₁`; the prime block intertwines only in the
two-isometry sense, with nonvanishing real characters. No `T*`-covariance. -/
def BT1Weak (D : BTSigmaDatum H Gu Gp) (Sd : ℝ → H ≃ₗᵢ[ℂ] H) (h : ℝ → ℝ) : Prop :=
  PCov D Sd ∧
  (∃ (Uu : ℝ → Gu →ₗᵢ[ℂ] Gu) (α₀ α₁ : ℝ → ℝ), ∀ W : ℝ, 0 ≤ W →
    D.Γ₀u (D.testVec h W) = (α₀ W : ℂ) • Uu W (D.Γ₀u (D.testVec h 0)) ∧
    D.Γ₁u (D.testVec h W) = (α₁ W : ℂ) • Uu W (D.Γ₁u (D.testVec h 0))) ∧
  (∃ (V U : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ),
    TwoIsometryPrime D h V U χ₀ χ₁ ∧ ∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0)

/-- `T*`-covariance along the orbit (the BT1 clause "`T` covariant"). -/
def TstarCov (D : BTSigmaDatum H Gu Gp) (Sd : ℝ → H ≃ₗᵢ[ℂ] H) (h : ℝ → ℝ) : Prop :=
  ∀ W b : ℝ, D.Tstar (D.testVec h (W + b)) = Sd b (D.Tstar (D.testVec h W))

/-- The unit clause of `BT1Weak` scales the unit block by `α₀α₁`. PROVED. -/
theorem bt1Weak_unitBlock (D : BTSigmaDatum H Gu Gp) (Sd : ℝ → H ≃ₗᵢ[ℂ] H) (h : ℝ → ℝ)
    (hB : BT1Weak D Sd h) :
    ∃ α : ℝ → ℝ, ∀ W : ℝ, 0 ≤ W → unitBlock D h W = (α W : ℂ) * unitBlock D h 0 := by
  obtain ⟨_, ⟨Uu, α₀, α₁, hu⟩, _⟩ := hB
  refine ⟨fun W => α₀ W * α₁ W, fun W hW => ?_⟩
  obtain ⟨e0, e1⟩ := hu W hW
  unfold unitBlock blockForm
  rw [e0, e1]
  simp only [inner_smul_left, inner_smul_right, LinearIsometry.inner_map_map, Complex.conj_ofReal]
  push_cast; ring

/-- **Growth survives the weakening.** PROVED (CovGrowth, verbatim). -/
theorem bt1Weak_growth (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) (hU : PUnit D) (Sd : ℝ → H ≃ₗᵢ[ℂ] H)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hB : BT1Weak D Sd (autocorr f))
    (R : ℝ) (hR : ∀ v, R < v → autocorr f v = 0) (W : ℝ) (hW0 : 0 ≤ W)
    (hW : R + Real.log N < W) :
    Real.exp (W / 2) * |∫ v : ℝ, autocorr f v * Real.exp (-v / 2)| ≤
      2 * (‖D.Tstar (D.testVec (autocorr f) W)‖ *
        ‖((D.testVec (autocorr f) 0 : D.Tstar.domain) : H)‖) :=
  pcov_graph_growth D S N hP hU Sd hB.1 f hf hc R hR W hW0 hW

/-- **Type `< 1/2` stays excluded.** PROVED (CovGrowth, verbatim). -/
theorem bt1Weak_type_lt_half (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) (hU : PUnit D) (Sd : ℝ → H ≃ₗᵢ[ℂ] H)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hB : BT1Weak D Sd (autocorr f))
    (R : ℝ) (hR : ∀ v, R < v → autocorr f v = 0) (M a : ℝ) (ha : a < 1 / 2)
    (hM : ∀ W : ℝ, 0 ≤ W → ‖D.Tstar (D.testVec (autocorr f) W)‖ ≤ M * Real.exp (a * W)) :
    ∫ v : ℝ, autocorr f v * Real.exp (-v / 2) = 0 :=
  pcov_type_lt_half_forces_zero D S N hP hU Sd hB.1 f hf hc R hR M a ha hM

/-- **Half-type results survive the weakening.** PROVED (CovHalfType, verbatim). -/
theorem bt1Weak_halfType (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) (hU : PUnit D) (Sd : ℝ → H ≃ₗᵢ[ℂ] H)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hB : BT1Weak D Sd (autocorr f))
    (R : ℝ) (hR : ∀ v, R < v → autocorr f v = 0) (M : ℝ)
    (hM : ∀ W : ℝ, 0 ≤ W → ‖D.Tstar (D.testVec (autocorr f) W)‖ ≤ M * Real.exp (W / 2))
    (W : ℝ) (hW0 : 0 ≤ W) (hW : R + Real.log N < W) :
    (Real.exp (W / 2) * |∫ v : ℝ, autocorr f v * Real.exp (-v / 2)| ≤
        2 * ‖((D.testVec (autocorr f) 0 : D.Tstar.domain) : H)‖ *
          ‖D.Tstar (D.testVec (autocorr f) W)‖ ∧
      ‖D.Tstar (D.testVec (autocorr f) W)‖ ≤ M * Real.exp (W / 2)) ∧
    (HalfHomogeneousEdge D (autocorr f) → ∀ W : ℝ, 0 ≤ W → stratumSN S N (autocorr f) W = 0) :=
  ⟨halfType_two_sided D S N hP hU Sd hB.1 f hf hc R hR M hM W hW0 hW,
    halfHomogeneous_forces_stratum_zero D S N hP hU f hf hc R hR⟩

/-- `T*`-covariance plus the (P-COV) orbit law make the edge motion Green-invariant. PROVED. -/
theorem tstarCov_greenInvariant (D : BTSigmaDatum H Gu Gp) (Sd : ℝ → H ≃ₗᵢ[ℂ] H)
    (hC : PCov D Sd) (h : ℝ → ℝ) (hT : TstarCov D Sd h) : GreenInvariantEdge D h := by
  intro W b _ _
  have hv := hC.2 h W b
  unfold greenForm
  rw [hT W b, hv]
  simp only [LinearIsometryEquiv.inner_map_map]

/-- **The obstruction still fires if `T*` stays covariant.** PROVED. -/
theorem bt1Weak_obstruction (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) (hU : PUnit D) (Sd : ℝ → H ≃ₗᵢ[ℂ] H)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hB : BT1Weak D Sd (autocorr f)) (hT : TstarCov D Sd (autocorr f))
    (R : ℝ) (hR : ∀ v, R < v → autocorr f v = 0) :
    ∫ v : ℝ, autocorr f v * Real.exp (-v / 2) = 0 :=
  greenInvariant_forces_unit_zero D S N hP hU f hf hc R hR
    (tstarCov_greenInvariant D Sd hB.1 _ hT)

end Datum

/-! ## (a) The two-isometry freedom admits the bump: an explicit `ℂ²` witness -/

/-- The prime boundary space of the witness. -/
abbrev E2 := EuclideanSpace ℂ (Fin 2)

/-- First basis vector. -/
noncomputable def e0 : E2 := EuclideanSpace.single 0 1
/-- Second basis vector. -/
noncomputable def e1 : E2 := EuclideanSpace.single 1 1

/-- `diag(l, 1)` on `ℂ²`, an isometry for `‖l‖ = 1`. -/
noncomputable def diagIso (l : ℂ) (hl : ‖l‖ = 1) : E2 →ₗᵢ[ℂ] E2 where
  toFun x := WithLp.toLp 2 (fun i => if i = 0 then l * x i else x i)
  map_add' x y := by
    ext i; by_cases hi : i = 0 <;> simp [hi, mul_add]
  map_smul' c x := by
    ext i; by_cases hi : i = 0 <;> simp [hi]; ring
  norm_map' x := by
    simp only [LinearMap.coe_mk, AddHom.coe_mk, EuclideanSpace.norm_eq, Fin.sum_univ_two]
    simp [hl]

theorem diagIso_e0 (l : ℂ) (hl : ‖l‖ = 1) : diagIso l hl e0 = l • e0 := by
  ext i; fin_cases i <;> simp [diagIso, e0]

theorem diagIso_e1 (l : ℂ) (hl : ‖l‖ = 1) : diagIso l hl e1 = e1 := by
  ext i; fin_cases i <;> simp [diagIso, e1]

theorem diagIso_one_apply (l : ℂ) (hl : ‖l‖ = 1) (h1 : l = 1) (x : E2) : diagIso l hl x = x := by
  subst h1; ext i; fin_cases i <;> simp [diagIso]

/-- `diag(l, 1)` is a scalar operator only if `l = 1`. PROVED. -/
theorem diagIso_nonscalar (l : ℂ) (hl : ‖l‖ = 1) (hl1 : l ≠ 1) :
    ¬ ∃ μ : ℂ, ∀ x : E2, diagIso l hl x = μ • x := by
  rintro ⟨μ, hμ⟩
  have h1 := hμ e1
  rw [diagIso_e1] at h1
  have hμ1 : μ = 1 := by
    have := congrArg (fun v : E2 => v 1) h1
    simp [e1] at this
    exact this.symm
  have h0 := hμ e0
  rw [diagIso_e0, hμ1, one_smul] at h0
  have := congrArg (fun v : E2 => v 0) h0
  simp [e0] at this
  exact hl1 this

/-- Witness phase parameter: `t(s) = −s/(2(1+|s|))`, `|t| < 1/2`. -/
noncomputable def tpar (s : ℝ) : ℝ := -s / (2 * (1 + |s|))

theorem tpar_sq_lt (s : ℝ) : tpar s ^ 2 < 1 := by
  have h1 : 0 < 1 + |s| := by positivity
  have h2 : |tpar s| < 1 := by
    unfold tpar
    rw [abs_div, abs_neg, abs_of_pos (by positivity : (0:ℝ) < 2 * (1 + |s|)),
      div_lt_one (by positivity)]
    linarith [abs_nonneg s]
  nlinarith [abs_nonneg (tpar s), sq_abs (tpar s)]

/-- Witness phase `λ(s) = √(1 − t²) + i t`, a unit complex number. -/
noncomputable def lam (s : ℝ) : ℂ := ⟨Real.sqrt (1 - tpar s ^ 2), tpar s⟩

theorem lam_norm (s : ℝ) : ‖lam s‖ = 1 := by
  have h := tpar_sq_lt s
  rw [Complex.norm_def, Complex.normSq_apply]
  simp only [lam]
  rw [← sq, ← sq, Real.sq_sqrt (by linarith)]
  simp

theorem lam_im (s : ℝ) : (lam s).im = tpar s := rfl

theorem lam_ne_one {s : ℝ} (hs : s ≠ 0) : lam s ≠ 1 := by
  intro h
  have := congrArg Complex.im h
  rw [lam_im, Complex.one_im] at this
  unfold tpar at this
  rw [div_eq_zero_iff] at this
  rcases this with h2 | h2
  · exact hs (by linarith)
  · have : 0 < 2 * (1 + |s|) := by positivity
    linarith

theorem relphase_norm (a b : ℝ) : ‖lam a * (starRingEnd ℂ) (lam b)‖ = 1 := by
  rw [norm_mul, Complex.norm_conj, lam_norm, lam_norm, one_mul]

theorem lam_mul_conj (s : ℝ) : lam s * (starRingEnd ℂ) (lam s) = 1 := by
  rw [Complex.mul_conj, Complex.normSq_eq_norm_sq, lam_norm]; simp

/-- The witness prime pair for a target value `s`: `Γ₀p = e₀`,
`Γ₁p = (1 + |s|) · diag(λ(s), 1) e₀`. -/
noncomputable def g1 (s : ℝ) : E2 := ((1 + |s| : ℝ) : ℂ) • diagIso (lam s) (lam_norm s) e0

/-- **Profile-level realisation.** The two-isometry boundary pair reproduces *any* real value
`s` of the prime block: `blockForm(e₀, g₁(s)) = I·s`. PROVED. -/
theorem witness_block (s : ℝ) : blockForm e0 (g1 s) e0 (g1 s) = I * (s : ℂ) := by
  unfold blockForm g1
  rw [diagIso_e0, smul_smul]
  simp only [inner_smul_left, inner_smul_right]
  have he : ⟪e0, e0⟫_ℂ = 1 := by simp [e0]
  rw [he, mul_one]
  apply Complex.ext
  · simp [lam]
  · simp only [Complex.sub_im, Complex.mul_im, Complex.conj_re, Complex.conj_im,
      Complex.ofReal_re, Complex.ofReal_im, Complex.I_re, Complex.I_im, map_mul,
      Complex.conj_ofReal]
    rw [lam_im]
    unfold tpar
    have h1 : (1 + |s|) ≠ 0 := by positivity
    field_simp
    simp only [Complex.one_re, Complex.one_im]
    ring

/-! ## Insertion-level re-test of a triplet -/

section Retest

variable {H Gu : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
  [NormedAddCommGroup Gu] [InnerProductSpace ℂ Gu]

/-- Boundary targets for the re-test: unit pair `(e, I·(T/2)·e)`, prime pair `(e₀, g₁(s))`. -/
noncomputable def retestTarget (S : Finset ℕ) (N : ℕ) (e : Gu) (h : ℝ → ℝ) (W : ℝ) :
    (Gu × E2) × (Gu × E2) :=
  ((e, e0), ((I * ((unitStratumTerm h W / 2 : ℝ) : ℂ)) • e, g1 (stratumSN S N h W)))

/-- The re-tested datum: same `T*`, same `Γ₀, Γ₁`, same Green identity and surjectivity; only
the test vectors are re-chosen (by surjectivity) to carry the targets. -/
noncomputable def retest (D : BTSigmaDatum H Gu E2) (S : Finset ℕ) (N : ℕ) (e : Gu) :
    BTSigmaDatum H Gu E2 :=
  { D with testVec := fun h W => Classical.choose (D.surj (retestTarget S N e h W)) }

theorem retest_spec (D : BTSigmaDatum H Gu E2) (S : Finset ℕ) (N : ℕ) (e : Gu) (h : ℝ → ℝ)
    (W : ℝ) :
    ((retest D S N e).Γ₀u ((retest D S N e).testVec h W) = e ∧
      (retest D S N e).Γ₀p ((retest D S N e).testVec h W) = e0) ∧
    ((retest D S N e).Γ₁u ((retest D S N e).testVec h W) =
        (I * ((unitStratumTerm h W / 2 : ℝ) : ℂ)) • e ∧
      (retest D S N e).Γ₁p ((retest D S N e).testVec h W) = g1 (stratumSN S N h W)) := by
  have hs := Classical.choose_spec (D.surj (retestTarget S N e h W))
  exact ⟨⟨congrArg (fun z => z.1.1) hs, congrArg (fun z => z.1.2) hs⟩,
    ⟨congrArg (fun z => z.2.1) hs, congrArg (fun z => z.2.2) hs⟩⟩

/-- The re-test satisfies finite-cut BT3 (for every test, every `W`). PROVED. -/
theorem retest_traceSN (D : BTSigmaDatum H Gu E2) (S : Finset ℕ) (N : ℕ) (e : Gu) :
    BTSigmaTraceSN (retest D S N e) S N := by
  intro f _ _ W _
  obtain ⟨⟨_, h0⟩, ⟨_, h1⟩⟩ := retest_spec D S N e (autocorr f) W
  unfold primeBlock
  rw [h0, h1, witness_block]

/-- The re-test satisfies (P-UNIT) when `e` is a unit vector. PROVED. -/
theorem retest_punit (D : BTSigmaDatum H Gu E2) (S : Finset ℕ) (N : ℕ) (e : Gu) (he : ‖e‖ = 1) :
    PUnit (retest D S N e) := by
  intro f _ _ W _
  obtain ⟨⟨u0, _⟩, ⟨u1, _⟩⟩ := retest_spec D S N e (autocorr f) W
  unfold unitBlock
  rw [u0, u1]
  unfold blockForm
  simp only [inner_smul_left, inner_smul_right, inner_self_eq_norm_sq_to_K, he]
  apply Complex.ext <;> (simp [map_ofNat]; try ring)

/-- Relative rotation of the witness along the orbit of `h`. -/
noncomputable def relRot (S : Finset ℕ) (N : ℕ) (h : ℝ → ℝ) (W : ℝ) : E2 →ₗᵢ[ℂ] E2 :=
  diagIso (lam (stratumSN S N h W) * (starRingEnd ℂ) (lam (stratumSN S N h 0)))
    (relphase_norm _ _)

/-- Real character on `Γ₁p`. -/
noncomputable def chi1 (S : Finset ℕ) (N : ℕ) (h : ℝ → ℝ) (W : ℝ) : ℝ :=
  (1 + |stratumSN S N h W|) / (1 + |stratumSN S N h 0|)

theorem chi1_pos (S : Finset ℕ) (N : ℕ) (h : ℝ → ℝ) (W : ℝ) : 0 < chi1 S N h W := by
  unfold chi1; positivity

/-- **Two-isometry covariance holds on the re-test, for every test `h`**, with `V_W = id`,
`U_W = diag(λ_W λ̄_0, 1)`, `χ₀ = 1`, `χ₁ > 0`. PROVED. -/
theorem retest_twoIsometry (D : BTSigmaDatum H Gu E2) (S : Finset ℕ) (N : ℕ) (e : Gu)
    (h : ℝ → ℝ) :
    TwoIsometryPrime (retest D S N e) h (fun _ => LinearIsometry.id) (relRot S N h)
      (fun _ => 1) (chi1 S N h) := by
  intro W _
  obtain ⟨⟨_, p0W⟩, ⟨_, p1W⟩⟩ := retest_spec D S N e h W
  obtain ⟨⟨_, p00⟩, ⟨_, p10⟩⟩ := retest_spec D S N e h 0
  refine ⟨?_, ?_⟩
  · rw [p0W, p00]; simp
  · rw [p1W, p10]
    unfold g1 relRot chi1
    simp only [diagIso_e0, LinearIsometry.map_smul, smul_smul]
    congr 1
    have hc := lam_mul_conj (stratumSN S N h 0)
    have hd : (1 + |stratumSN S N h 0|) ≠ 0 := by positivity
    generalize lam (stratumSN S N h W) = a at *
    generalize lam (stratumSN S N h 0) = b at *
    generalize (1 + |stratumSN S N h W|) = A at *
    generalize (1 + |stratumSN S N h 0|) = B at *
    have hB : (B : ℂ) ≠ 0 := by exact_mod_cast hd
    rw [Complex.ofReal_div]
    field_simp
    linear_combination (-(A : ℂ) * a) * hc

/-- **Separation.** On the re-test of any triplet with `Gp = ℂ²` and a unit vector in `Gu`,
with a prime `p ∈ S`, `p ≤ N`: finite-cut BT3 and (P-UNIT) hold; same-isometry covariance fails
for some smooth compactly supported test (BoundaryPW); two-isometry covariance holds on that
test; and there the relative operator `V_W* U_W = U_W` is non-scalar at some `W ≥ 0`.
PROVED (insertion level; (P-COV) not claimed). -/
theorem retest_separates (D : BTSigmaDatum H Gu E2) (S : Finset ℕ) (N : ℕ) (e : Gu)
    (he : ‖e‖ = 1) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S) (hpN : p ≤ N) :
    BTSigmaTraceSN (retest D S N e) S N ∧ PUnit (retest D S N e) ∧
    ∃ f : ℝ → ℝ, ContDiff ℝ ∞ f ∧ HasCompactSupport f ∧
      (∀ (U : ℝ → E2 →ₗᵢ[ℂ] E2) (χ₀ χ₁ : ℝ → ℝ), (∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0) →
        ¬ OrbitCovariantPrime (retest D S N e) (autocorr f) U χ₀ χ₁) ∧
      TwoIsometryPrime (retest D S N e) (autocorr f) (fun _ => LinearIsometry.id)
        (relRot S N (autocorr f)) (fun _ => 1) (chi1 S N (autocorr f)) ∧
      ∃ W : ℝ, 0 ≤ W ∧ ¬ ∃ μ : ℂ, ∀ x : E2, relRot S N (autocorr f) W x = μ • x := by
  have hP := retest_traceSN D S N e
  refine ⟨hP, retest_punit D S N e he, ?_⟩
  obtain ⟨f, hf, hc, hf0, hf1⟩ := exists_smooth_test
  obtain ⟨R, hR⟩ := autocorr_eventually_zero hc
  have hpos := stratumSN_pos S N hp hpS hpN hf.continuous hc hf0 hf1
  refine ⟨f, hf, hc, ?_, retest_twoIsometry D S N e _, ?_⟩
  · intro U χ₀ χ₁ hχ hcov
    have hz := covariant_prime_forces_stratum_zero (retest D S N e) S N hP f hf hc R hR U χ₀ χ₁
      hcov hχ (Real.log p) (Real.log_natCast_nonneg p)
    linarith
  · by_contra hall
    push Not at hall
    have hid : ∀ W : ℝ, 0 ≤ W → ∀ x, relRot S N (autocorr f) W x = x := by
      intro W hW x
      obtain ⟨μ, hμ⟩ := hall W hW
      have hμ1 : lam (stratumSN S N (autocorr f) W) *
          (starRingEnd ℂ) (lam (stratumSN S N (autocorr f) 0)) = 1 := by
        by_contra hl1
        exact diagIso_nonscalar _ (relphase_norm _ _) hl1 ⟨μ, hμ⟩
      exact diagIso_one_apply _ _ hμ1 x
    have hrel : RelPhaseTrivial (retest D S N e) (autocorr f) (fun _ => LinearIsometry.id)
        (relRot S N (autocorr f)) := by
      intro W hW
      simp only [LinearIsometry.id_apply]
      rw [hid W hW]
    have hz := twoIsometry_relPhase_forces_stratum_zero (retest D S N e) S N hP f hf hc R hR
      _ _ _ _ (retest_twoIsometry D S N e _) hrel
      (fun W _ => by have := chi1_pos S N (autocorr f) W; simp only [one_mul]; exact this.ne')
      (Real.log p) (Real.log_natCast_nonneg p)
    linarith

end Retest

end GppBTSigmaTwoIsometry

#print axioms GppBTSigmaTwoIsometry.covariant_iff_twoIsometry_diag
#print axioms GppBTSigmaTwoIsometry.primeBlock_twoIsometry
#print axioms GppBTSigmaTwoIsometry.scaling_forces_stratum_zero
#print axioms GppBTSigmaTwoIsometry.twoIsometry_relPhase_forces_stratum_zero
#print axioms GppBTSigmaTwoIsometry.relPhase_of_same
#print axioms GppBTSigmaTwoIsometry.bt1Weak_unitBlock
#print axioms GppBTSigmaTwoIsometry.bt1Weak_growth
#print axioms GppBTSigmaTwoIsometry.bt1Weak_type_lt_half
#print axioms GppBTSigmaTwoIsometry.bt1Weak_halfType
#print axioms GppBTSigmaTwoIsometry.tstarCov_greenInvariant
#print axioms GppBTSigmaTwoIsometry.bt1Weak_obstruction
#print axioms GppBTSigmaTwoIsometry.diagIso_e0
#print axioms GppBTSigmaTwoIsometry.diagIso_e1
#print axioms GppBTSigmaTwoIsometry.diagIso_one_apply
#print axioms GppBTSigmaTwoIsometry.diagIso_nonscalar
#print axioms GppBTSigmaTwoIsometry.tpar_sq_lt
#print axioms GppBTSigmaTwoIsometry.lam_norm
#print axioms GppBTSigmaTwoIsometry.lam_im
#print axioms GppBTSigmaTwoIsometry.lam_ne_one
#print axioms GppBTSigmaTwoIsometry.relphase_norm
#print axioms GppBTSigmaTwoIsometry.lam_mul_conj
#print axioms GppBTSigmaTwoIsometry.witness_block
#print axioms GppBTSigmaTwoIsometry.retest_spec
#print axioms GppBTSigmaTwoIsometry.retest_traceSN
#print axioms GppBTSigmaTwoIsometry.retest_punit
#print axioms GppBTSigmaTwoIsometry.chi1_pos
#print axioms GppBTSigmaTwoIsometry.retest_twoIsometry
#print axioms GppBTSigmaTwoIsometry.retest_separates
