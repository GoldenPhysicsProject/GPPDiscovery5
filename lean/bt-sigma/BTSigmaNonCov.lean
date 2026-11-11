import GppVerify.RiemannHypothesis.BTSigmaBoundaryPW
import GppVerify.RiemannHypothesis.BTSigmaCovHalfType

/-!
# BT-Σ: the (P-NONCOV) test

Scratch module (Heart of the Cosmos / Grok, GPPDiscovery5 seat) for Turbo's question on
GPP-bridge PR #11, comment 6804139713.

**Part 1 (scope gap of `BTSigmaBoundaryPW`).** Independent isometries `U_W` on `Γ₀p` and `V_W` on
`Γ₁p`. `primeBlock_indep`: the prime block becomes `χ₀χ₁ · blockForm(U_W a, V_W b)`, which depends
on the relative pairing `⟪U_W a, V_W b⟫`. `indep_kill_of_pairing`: the kill survives when that
pairing is preserved (`V_W = U_W` is the old case; the exact condition is
`Im⟪U_W a, V_W b⟫ = k(W)·Im⟪a, b⟫` with real `k`, a sufficient form of which is proved). `phase_evader`: otherwise it fails. On
`Gp = ℂ` with `U_W = id`, `χ ≡ 1` and `V_W` the phase `e^{i(θ(W) − θ(0))}`, where
`θ = arcsin(−g/(2K))`, the prime block equals `I·g(W)` for **any** real profile with `|g| ≤ 2K`.
That includes a compactly supported bump. Even a scalar relative phase evades the kill, provided
it is not real.

**Part 2 (P-NONCOV).** `PNonCov` is a named Prop. `pnoncov_of_pcov_bt3`: under BT3 the
non-covariance clause is a theorem, not a hypothesis. So `P-NONCOV ⇔ P-COV ∧ BT3` on every
nonnegative nonzero smooth test, and (b) fails in the strongest form: BT3 forces the prime block
**not** to intertwine (single-isometry sense). `bounded_unit_obstruction`, `finiteDim_obstruction`:
P-COV + P-UNIT + boundedness of `Γ₀u, Γ₁u` in the `H`-norm is contradictory for any test with a
nonzero unit integral. In particular no finite-dimensional toy exists, so the toy asked for in (a)
cannot be finite. A (P-NONCOV) datum needs `Γ_u` that is unbounded in the `H`-norm (graph-bounded
only) on an infinite-dimensional `H`. Its existence is the existing OPEN S_b problem.

All PROVED, standard axioms. No `sorry`, no `admit`, no custom axiom. S_b existence, (P-COV),
(P-UNIT), (P-LOG), (BT-Σ), ARP and RH remain OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo.
-/

namespace GppBTSigmaNonCov

open GppBTSigmaLocalEuler GppBTSigma GppArithmeticReflection Complex MeasureTheory
open GppBTSigmaCovGrowth GppBTSigmaCovObstruction GppBTSigmaCovHalfType GppBTSigmaBoundaryPW
open scoped ContDiff

section Datum

variable {H Gu Gp : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
  [NormedAddCommGroup Gu] [InnerProductSpace ℂ Gu]
  [NormedAddCommGroup Gp] [InnerProductSpace ℂ Gp]

/-- Covariance along the orbit with **independent** isometries `U_W` on `Γ₀p` and `V_W` on `Γ₁p`. -/
def OrbitIndepCovariantPrime (D : BTSigmaDatum H Gu Gp) (h : ℝ → ℝ) (U V : ℝ → Gp →ₗᵢ[ℂ] Gp)
    (χ₀ χ₁ : ℝ → ℝ) : Prop :=
  ∀ W : ℝ, 0 ≤ W →
    D.Γ₀p (D.testVec h W) = (χ₀ W : ℂ) • U W (D.Γ₀p (D.testVec h 0)) ∧
    D.Γ₁p (D.testVec h W) = (χ₁ W : ℂ) • V W (D.Γ₁p (D.testVec h 0))

/-- Under independent isometries the prime block is `χ₀χ₁ · blockForm(U_W a, V_W b)`. PROVED. -/
theorem primeBlock_indep (D : BTSigmaDatum H Gu Gp) (h : ℝ → ℝ) (U V : ℝ → Gp →ₗᵢ[ℂ] Gp)
    (χ₀ χ₁ : ℝ → ℝ) (hcov : OrbitIndepCovariantPrime D h U V χ₀ χ₁) (W : ℝ) (hW : 0 ≤ W) :
    primeBlock D h W = ((χ₀ W * χ₁ W : ℝ) : ℂ) *
      blockForm (U W (D.Γ₀p (D.testVec h 0))) (V W (D.Γ₁p (D.testVec h 0)))
        (U W (D.Γ₀p (D.testVec h 0))) (V W (D.Γ₁p (D.testVec h 0))) := by
  obtain ⟨e0, e1⟩ := hcov W hW
  unfold primeBlock blockForm
  rw [e0, e1]
  simp only [inner_smul_left, inner_smul_right, Complex.conj_ofReal]
  push_cast; ring

/-- A prime block that scales by a nonvanishing real factor forces `stratumSN ≡ 0` under BT3. PROVED. -/
theorem scaled_primeBlock_forces_stratum_zero (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (R : ℝ) (hR : ∀ v, R < v → autocorr f v = 0) (k : ℝ → ℝ) (hk : ∀ W : ℝ, 0 ≤ W → k W ≠ 0)
    (hs : ∀ W : ℝ, 0 ≤ W → primeBlock D (autocorr f) W = (k W : ℂ) * primeBlock D (autocorr f) 0) :
    ∀ W : ℝ, 0 ≤ W → stratumSN S N (autocorr f) W = 0 := by
  have key : ∀ W : ℝ, 0 ≤ W →
      stratumSN S N (autocorr f) W = k W * stratumSN S N (autocorr f) 0 := by
    intro W hW
    have e := hs W hW
    rw [hP f hf hc W hW, hP f hf hc 0 le_rfl] at e
    have e2 : ((stratumSN S N (autocorr f) W : ℝ) : ℂ) =
        ((k W * stratumSN S N (autocorr f) 0 : ℝ) : ℂ) := by
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
    · exact absurd h1 (hk W₀ hW₀)
    · exact h1
  intro W hW
  rw [key W hW, h0, mul_zero]

/-- **Where the kill survives.** If the relative pairing `⟪U_W a, V_W b⟫ = ⟪a, b⟫` is preserved
along the orbit (e.g. `V_W = U_W`), independent covariance still forces `stratumSN ≡ 0`. PROVED. -/
theorem indep_kill_of_pairing (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (R : ℝ) (hR : ∀ v, R < v → autocorr f v = 0)
    (U V : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ)
    (hcov : OrbitIndepCovariantPrime D (autocorr f) U V χ₀ χ₁)
    (hχ : ∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0)
    (hpair : ∀ W : ℝ, 0 ≤ W →
      inner ℂ (U W (D.Γ₀p (D.testVec (autocorr f) 0))) (V W (D.Γ₁p (D.testVec (autocorr f) 0))) =
        inner ℂ (D.Γ₀p (D.testVec (autocorr f) 0)) (D.Γ₁p (D.testVec (autocorr f) 0))) :
    ∀ W : ℝ, 0 ≤ W → stratumSN S N (autocorr f) W = 0 := by
  refine scaled_primeBlock_forces_stratum_zero D S N hP f hf hc R hR (fun W => χ₀ W * χ₁ W) hχ ?_
  intro W hW
  rw [primeBlock_indep D _ U V χ₀ χ₁ hcov W hW]
  congr 1
  unfold primeBlock blockForm
  rw [← inner_conj_symm (V W _) (U W _), hpair W hW, inner_conj_symm]

end Datum

/-! ## The evading intertwiner: a non-real relative phase on `Gp = ℂ` -/

/-- Multiplication by the phase `e^{iθ}` as a `ℂ`-linear isometry of `ℂ`. -/
noncomputable def phaseIso (θ : ℝ) : ℂ →ₗᵢ[ℂ] ℂ where
  toLinearMap := LinearMap.lsmul ℂ ℂ (cexp (θ * I))
  norm_map' x := by
    simp [Complex.norm_exp_ofReal_mul_I]

theorem phaseIso_apply (θ : ℝ) (z : ℂ) : phaseIso θ z = cexp (θ * I) * z := rfl

/-- The zero phase is the identity. PROVED. -/
theorem phaseIso_zero (z : ℂ) : phaseIso 0 z = z := by
  simp [phaseIso_apply]

/-- `blockForm(1, K e^{iθ}) = I · (−2K sin θ)`. PROVED. -/
theorem phase_blockForm (K θ : ℝ) :
    blockForm (1 : ℂ) ((K : ℂ) * cexp (θ * I)) 1 ((K : ℂ) * cexp (θ * I)) =
      I * ((-2 * K * Real.sin θ : ℝ) : ℂ) := by
  unfold blockForm
  apply Complex.ext <;>
    simp [Complex.exp_ofReal_mul_I_re, Complex.exp_ofReal_mul_I_im, Complex.sin_ofReal_re] <;> ring

/-- **The kill fails for independent isometries.** With `U_W = id`, `χ₀ = χ₁ = 1`, boundary values
`a = 1`, `b = K e^{iθ(0)}` at `W = 0` and `V_W = phaseIso(θ(W) − θ(0))` (so `V_0 = id`), where
`θ = arcsin(−g/(2K))`, the prime block along the orbit is `I·g(W)` for any real profile `g` with
`|g| ≤ 2K`, for example a compactly supported bump such as `stratumSN`. PROVED. -/
theorem phase_evader (g : ℝ → ℝ) (K : ℝ) (hK : 0 < K) (hg : ∀ W, |g W| ≤ 2 * K) (W : ℝ) :
    let θ : ℝ → ℝ := fun W => Real.arcsin (-g W / (2 * K))
    let b : ℂ := (K : ℂ) * cexp (θ 0 * I)
    blockForm (1 : ℂ) (phaseIso (θ W - θ 0) b) 1 (phaseIso (θ W - θ 0) b) = I * (g W : ℂ) := by
  intro θ b
  have hb : phaseIso (θ W - θ 0) b = (K : ℂ) * cexp (θ W * I) := by
    simp only [phaseIso_apply, b]
    rw [mul_left_comm, ← Complex.exp_add]
    congr 2; push_cast; ring
  rw [hb, phase_blockForm]
  congr 1
  have h1 : -1 ≤ -g W / (2 * K) := by
    rw [le_div_iff₀ (by linarith)]; have := (abs_le.mp (hg W)).2; linarith
  have h2 : -g W / (2 * K) ≤ 1 := by
    rw [div_le_iff₀ (by linarith)]; have := (abs_le.mp (hg W)).1; linarith
  have hs : Real.sin (θ W) = -g W / (2 * K) := Real.sin_arcsin h1 h2
  rw [hs]
  congr 1
  field_simp

/-! ## Part 2: (P-NONCOV) -/

section Datum2

variable {H Gu Gp : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
  [NormedAddCommGroup Gu] [InnerProductSpace ℂ Gu]
  [NormedAddCommGroup Gp] [InnerProductSpace ℂ Gp]

/-- **(P-NONCOV), named, not asserted.** On the fixed test orbit of `f`: a dilation `S_b`
(P-COV) exists; the prime block is the prescribed profile `I·stratumSN(W)` exactly; and `Γ₀p, Γ₁p`
see `W` non-covariantly. That means there is no isometry `U_W` and no real nonvanishing character
with `Γ_ip ψ_W = χ_i(W)·U_W Γ_ip ψ_0`. -/
def PNonCov (D : BTSigmaDatum H Gu Gp) (Sd : ℝ → H ≃ₗᵢ[ℂ] H) (S : Finset ℕ) (N : ℕ)
    (f : ℝ → ℝ) : Prop :=
  PCov D Sd ∧
  (∀ W : ℝ, 0 ≤ W → primeBlock D (autocorr f) W = I * (stratumSN S N (autocorr f) W : ℂ)) ∧
  (∀ (U : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ), (∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0) →
    ¬ OrbitCovariantPrime D (autocorr f) U χ₀ χ₁)

/-- **The non-covariance clause is automatic.** For any smooth, compactly supported, nonnegative,
nonzero test, P-COV and finite-cut BT3 (with a prime `p ∈ S`, `p ≤ N`) already give
(P-NONCOV). So BT3 forces the prime block **not** to intertwine, the opposite of (b). PROVED. -/
theorem pnoncov_of_pcov_bt3 (D : BTSigmaDatum H Gu Gp) (Sd : ℝ → H ≃ₗᵢ[ℂ] H) (S : Finset ℕ)
    (N : ℕ) (hC : PCov D Sd) (hP : BTSigmaTraceSN D S N) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S)
    (hpN : p ≤ N) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hf0 : ∀ x, 0 ≤ f x) {x₀ : ℝ} (hx : f x₀ ≠ 0) :
    PNonCov D Sd S N f := by
  obtain ⟨R, hR⟩ := autocorr_eventually_zero hc
  refine ⟨hC, fun W hW => hP f hf hc W hW, fun U χ₀ χ₁ hχ hcov => ?_⟩
  have hz := covariant_prime_forces_stratum_zero D S N hP f hf hc R hR U χ₀ χ₁ hcov hχ
    (Real.log p) (Real.log_natCast_nonneg p)
  have hpos := stratumSN_pos S N hp hpS hpN hf.continuous hc hf0 hx
  linarith

/-- (P-NONCOV) contains P-COV and BT3 on the orbit, so `P-NONCOV ⇔ P-COV ∧ BT3(f)` for such
tests. PROVED. -/
theorem pnoncov_iff (D : BTSigmaDatum H Gu Gp) (Sd : ℝ → H ≃ₗᵢ[ℂ] H) (S : Finset ℕ)
    (N : ℕ) (hP : BTSigmaTraceSN D S N) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S)
    (hpN : p ≤ N) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hf0 : ∀ x, 0 ≤ f x) {x₀ : ℝ} (hx : f x₀ ≠ 0) :
    PNonCov D Sd S N f ↔ PCov D Sd :=
  ⟨fun h => h.1, fun hC => pnoncov_of_pcov_bt3 D Sd S N hC hP hp hpS hpN f hf hc hf0 hx⟩

/-- **Bounded unit boundary maps are impossible.** P-COV + P-UNIT + `‖Γ₀u x‖ ≤ C₀‖x‖`,
`‖Γ₁u x‖ ≤ C₁‖x‖` contradict a nonzero unit integral `∫ (f⋆f̃)(v) e^{-v/2} dv`. PROVED. -/
theorem bounded_unit_obstruction (D : BTSigmaDatum H Gu Gp) (Sd : ℝ → H ≃ₗᵢ[ℂ] H)
    (hC : PCov D Sd) (hU : PUnit D) (C₀ C₁ : ℝ)
    (h0 : ∀ x : D.Tstar.domain, ‖D.Γ₀u x‖ ≤ C₀ * ‖(x : H)‖)
    (h1 : ∀ x : D.Tstar.domain, ‖D.Γ₁u x‖ ≤ C₁ * ‖(x : H)‖)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (hI : (∫ v : ℝ, autocorr f v * Real.exp (-v / 2)) ≠ 0) : False := by
  set c := |∫ v : ℝ, autocorr f v * Real.exp (-v / 2)| with hcdef
  have hcpos : 0 < c := abs_pos.mpr hI
  set n := ‖((D.testVec (autocorr f) 0 : D.Tstar.domain) : H)‖
  set M := 2 * ((max C₀ 0) * n * ((max C₁ 0) * n)) with hM
  have hMnn : 0 ≤ M := by positivity
  set W : ℝ := 2 * (M / c) with hWdef
  have hW : 0 ≤ W := by positivity
  have lower := punit_unit_product_lower D hU f hf hc W hW
  have hn : ‖((D.testVec (autocorr f) W : D.Tstar.domain) : H)‖ = n :=
    pcov_norm_const D Sd hC _ W
  have b0 : ‖D.Γ₀u (D.testVec (autocorr f) W)‖ ≤ max C₀ 0 * n := by
    have := h0 (D.testVec (autocorr f) W); rw [hn] at this
    exact this.trans (mul_le_mul_of_nonneg_right (le_max_left _ _) (norm_nonneg _))
  have b1 : ‖D.Γ₁u (D.testVec (autocorr f) W)‖ ≤ max C₁ 0 * n := by
    have := h1 (D.testVec (autocorr f) W); rw [hn] at this
    exact this.trans (mul_le_mul_of_nonneg_right (le_max_left _ _) (norm_nonneg _))
  have upper : 2 * (‖D.Γ₀u (D.testVec (autocorr f) W)‖ * ‖D.Γ₁u (D.testVec (autocorr f) W)‖)
      ≤ M := by
    rw [hM]
    have := mul_le_mul b0 b1 (norm_nonneg _) (by positivity)
    linarith
  have hexp : W / 2 + 1 ≤ Real.exp (W / 2) := Real.add_one_le_exp _
  have : (M / c + 1) * c ≤ Real.exp (W / 2) * c := by
    apply mul_le_mul_of_nonneg_right _ hcpos.le
    rw [hWdef] at hexp; linarith
  have e : (M / c + 1) * c = M + c := by field_simp
  linarith

/-- **No finite-dimensional toy.** If `H` is finite-dimensional, P-COV + P-UNIT contradict a
nonzero unit integral. A (P-NONCOV) datum with P-UNIT needs `Γ_u` unbounded in the `H`-norm, on an
infinite-dimensional `H`. PROVED. -/
theorem finiteDim_obstruction [FiniteDimensional ℂ H]
    (D : BTSigmaDatum H Gu Gp) (Sd : ℝ → H ≃ₗᵢ[ℂ] H)
    (hC : PCov D Sd) (hU : PUnit D) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f)
    (hc : HasCompactSupport f) (hI : (∫ v : ℝ, autocorr f v * Real.exp (-v / 2)) ≠ 0) :
    False := by
  let L0 := LinearMap.toContinuousLinearMap D.Γ₀u
  let L1 := LinearMap.toContinuousLinearMap D.Γ₁u
  exact bounded_unit_obstruction D Sd hC hU ‖L0‖ ‖L1‖ (fun x => L0.le_opNorm x)
    (fun x => L1.le_opNorm x) f hf hc hI

end Datum2

end GppBTSigmaNonCov

#print axioms GppBTSigmaNonCov.primeBlock_indep
#print axioms GppBTSigmaNonCov.scaled_primeBlock_forces_stratum_zero
#print axioms GppBTSigmaNonCov.indep_kill_of_pairing
#print axioms GppBTSigmaNonCov.phaseIso_zero
#print axioms GppBTSigmaNonCov.phase_blockForm
#print axioms GppBTSigmaNonCov.phase_evader
#print axioms GppBTSigmaNonCov.pnoncov_of_pcov_bt3
#print axioms GppBTSigmaNonCov.pnoncov_iff
#print axioms GppBTSigmaNonCov.bounded_unit_obstruction
#print axioms GppBTSigmaNonCov.finiteDim_obstruction
