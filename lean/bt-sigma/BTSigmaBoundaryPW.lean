import GppVerify.RiemannHypothesis.BTSigmaCovPerturbation

/-!
# BT-Σ: the boundary-carried perturbation `P_W` — form-level carriage is automatic, covariant carriage is impossible

Scratch module (Heart of the Cosmos / Grok, GPPDiscovery5 seat) for Turbo's question on
GPP-bridge PR #11, comment 6803775354: construct, or prove impossible, a boundary potential at the
cuts `W = log p^k`, carried by the prime-block boundary maps `Γ₀p, Γ₁p`, of weight
`halfWeight(p^k)` and profile `h`, whose Green form on the test orbit is exactly `I·stratumSN(W)`
(the bulk toy's profile from `BTSigmaCovPerturbation`, rebuilt as boundary data).

**Scope.** Fixed test orbit `ψ_W = testVec (autocorr f) W`; no zeros used as input; nothing about
existence of `S_b` for a genuine `BTSigmaDatum` beyond what is proved below.

**Results (all PROVED, standard axioms).**
* `boundary_carriage_automatic`: under (P-COV), finite-cut BT3, (P-UNIT) and any orbit
  decomposition `T* ψ_W = S_W (e^{W/2} A ψ_0 + P_W ψ_0)` with late-vanishing perturbation form,
  `form_{P_W}(ψ_0) = primeBlock(W)` and `e^{W/2}·form_A(ψ_0) = unitBlock(W)` for all `W ≥ 0`.
  At the level of forms the perturbation **is** carried by `Γ₀p, Γ₁p` and the fixed part by
  `Γ₀u, Γ₁u`; "carried by the prime block" adds nothing beyond BT3 along the orbit.
* `punit_unitBlock_covariant`: the unit block scales by the character `e^{W/2}` on the orbit, so
  the unit stratum **is** compatible with a covariant (BT1-intertwining) boundary map.
* `primeBlock_covariant`, `covariant_prime_forces_stratum_zero`: if the prime boundary maps are
  covariant along the orbit, `Γ₀p ψ_W = χ₀(W)·U_W Γ₀p ψ_0`, `Γ₁p ψ_W = χ₁(W)·U_W Γ₁p ψ_0` with
  `U_W` a linear isometry of `Gp` and real `χ₀χ₁ ≠ 0` (the BT1 clause "Γ₀, Γ₁ intertwining",
  allowing any real character), then `primeBlock(W) = χ₀χ₁(W)·primeBlock(0)` and finite-cut BT3
  forces `stratumSN ≡ 0` on `W ≥ 0`.
* `stratumSN_pos`, `exists_smooth_test`, `autocorr_eventually_zero`: for `f ≥ 0` continuous,
  compactly supported, `f ≠ 0`, and a prime `p ∈ S` with `p ≤ N`,
  `stratumSN(log p) ≥ halfWeight(p)·‖f‖² > 0`; a smooth such `f` exists (a `ContDiffBump`).
* `no_covariant_prime_boundary`: hence for every datum with finite-cut BT3 and `S` containing a
  prime `p ≤ N` there is a smooth compactly supported test `f` for which **no** covariant prime
  boundary exists along its orbit. Unconditional.
* `translate_jumpset_empty`: the concrete local-trace model fails the same way. If the prime
  boundary values are point traces (jumps) at a finite fixed cut set `C` and `S_W` is a translation,
  the jump set `J` of `ψ_0` must satisfy `J + W ⊆ C` for all `W ≥ 0`, so `J = ∅`: every orbit
  vector is jump-free and the prime block vanishes.

**Answer: IMPOSSIBLE for a covariant boundary potential; otherwise it needs a further hypothesis
that contradicts BT1.** The step where bulk-to-boundary breaks is exactly the intertwining clause:
the toy's profile `∫ φ(x) f(x − W) dx` samples the bulk potential at every `x` while the density
moves; a boundary channel samples one fixed cut, and covariance turns its value into a character
`χ(W)`, which cannot be a compactly supported bump. A boundary carrier must therefore make
`Γ₀p, Γ₁p` see `W` non-covariantly, i.e. the prime block must not intertwine with `S_b` (while the
unit block may). Whether such a non-covariant genuine datum exists is **OPEN**; S_b existence,
(P-COV), (P-UNIT), (P-LOG), (P-θ), (P-LIM), (BT-Σ), ARP and RH remain **OPEN**.

No `sorry`, no `admit`, no custom axiom.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo.
-/

namespace GppBTSigmaBoundaryPW

open GppBTSigmaLocalEuler GppBTSigma GppArithmeticReflection Complex MeasureTheory
open GppBTSigmaCovGrowth GppBTSigmaCovObstruction GppBTSigmaCovPerturbation GppBTSigmaC1
open scoped ContDiff

section Datum

variable {H Gu Gp : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
  [NormedAddCommGroup Gu] [InnerProductSpace ℂ Gu]
  [NormedAddCommGroup Gp] [InnerProductSpace ℂ Gp]

/-- **Form-level carriage is automatic.** The forced perturbation form is the prime block and the
fixed part is the unit block, for every `W ≥ 0`. PROVED. -/
theorem boundary_carriage_automatic (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) (hU : PUnit D) (Sd : ℝ → H ≃ₗᵢ[ℂ] H) (hC : PCov D Sd)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (R : ℝ) (hR : ∀ v, R < v → autocorr f v = 0)
    (A : H →ₗ[ℂ] H) (P : ℝ → H →ₗ[ℂ] H)
    (hdec : ∀ W : ℝ, 0 ≤ W → D.Tstar (D.testVec (autocorr f) W) =
      Sd W ((Real.exp (W / 2) : ℂ) • A (D.testVec (autocorr f) 0 : H) +
        P W (D.testVec (autocorr f) 0 : H)))
    (W₁ : ℝ) (hlate : ∀ W : ℝ, W₁ ≤ W → opForm (P W) (D.testVec (autocorr f) 0 : H) = 0) :
    ∀ W : ℝ, 0 ≤ W →
      opForm (P W) (D.testVec (autocorr f) 0 : H) = primeBlock D (autocorr f) W ∧
      (Real.exp (W / 2) : ℂ) * opForm A (D.testVec (autocorr f) 0 : H) =
        unitBlock D (autocorr f) W := by
  intro W hW
  obtain ⟨hA, hPW⟩ := perturbation_forced D S N hP hU Sd hC f hf hc R hR A P hdec W₁ hlate
  refine ⟨?_, ?_⟩
  · rw [hPW W hW, hP f hf hc W hW]
  · rw [hA, hU f hf hc W hW, unitStratumTerm]
    push_cast; ring

/-- **The unit block is covariant with character `e^{W/2}`.** PROVED. -/
theorem punit_unitBlock_covariant (D : BTSigmaDatum H Gu Gp) (hU : PUnit D)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f) (W : ℝ) (hW : 0 ≤ W) :
    unitBlock D (autocorr f) W = (Real.exp (W / 2) : ℂ) * unitBlock D (autocorr f) 0 := by
  rw [hU f hf hc W hW, hU f hf hc 0 le_rfl, unitStratumTerm, unitStratumTerm, zero_div,
    Real.exp_zero, one_mul]
  push_cast; ring

/-- Covariant prime boundary along the orbit: the BT1 intertwining clause on the prime block,
with a linear isometry `U_W` of `Gp` and real characters `χ₀, χ₁`. -/
def OrbitCovariantPrime (D : BTSigmaDatum H Gu Gp) (h : ℝ → ℝ) (U : ℝ → Gp →ₗᵢ[ℂ] Gp)
    (χ₀ χ₁ : ℝ → ℝ) : Prop :=
  ∀ W : ℝ, 0 ≤ W →
    D.Γ₀p (D.testVec h W) = (χ₀ W : ℂ) • U W (D.Γ₀p (D.testVec h 0)) ∧
    D.Γ₁p (D.testVec h W) = (χ₁ W : ℂ) • U W (D.Γ₁p (D.testVec h 0))

/-- A covariant prime boundary scales the prime block by `χ₀χ₁`. PROVED. -/
theorem primeBlock_covariant (D : BTSigmaDatum H Gu Gp) (h : ℝ → ℝ) (U : ℝ → Gp →ₗᵢ[ℂ] Gp)
    (χ₀ χ₁ : ℝ → ℝ) (hcov : OrbitCovariantPrime D h U χ₀ χ₁) (W : ℝ) (hW : 0 ≤ W) :
    primeBlock D h W = ((χ₀ W * χ₁ W : ℝ) : ℂ) * primeBlock D h 0 := by
  obtain ⟨e0, e1⟩ := hcov W hW
  unfold primeBlock blockForm
  rw [e0, e1]
  simp only [inner_smul_left, inner_smul_right, LinearIsometry.inner_map_map, Complex.conj_ofReal]
  push_cast; ring

/-- **Covariant prime carriage kills the stratum.** Finite-cut BT3 + a covariant prime boundary
with nonvanishing character force `stratumSN ≡ 0` on `W ≥ 0`. PROVED. -/
theorem covariant_prime_forces_stratum_zero (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (R : ℝ) (hR : ∀ v, R < v → autocorr f v = 0)
    (U : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ)
    (hcov : OrbitCovariantPrime D (autocorr f) U χ₀ χ₁)
    (hχ : ∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0) :
    ∀ W : ℝ, 0 ≤ W → stratumSN S N (autocorr f) W = 0 := by
  have key : ∀ W : ℝ, 0 ≤ W →
      stratumSN S N (autocorr f) W = χ₀ W * χ₁ W * stratumSN S N (autocorr f) 0 := by
    intro W hW
    have e := primeBlock_covariant D _ U χ₀ χ₁ hcov W hW
    rw [hP f hf hc W hW, hP f hf hc 0 le_rfl] at e
    have e2 : ((stratumSN S N (autocorr f) W : ℝ) : ℂ) =
        ((χ₀ W * χ₁ W * stratumSN S N (autocorr f) 0 : ℝ) : ℂ) := by
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
    · exact absurd h1 (hχ W₀ hW₀)
    · exact h1
  intro W hW
  rw [key W hW, h0, mul_zero]

end Datum

/-! ## The stratum is not identically zero -/

theorem autocorr_nonneg {f : ℝ → ℝ} (hf0 : ∀ x, 0 ≤ f x) (v : ℝ) : 0 ≤ autocorr f v := by
  unfold autocorr
  exact integral_nonneg fun y => mul_nonneg (hf0 _) (hf0 _)

theorem autocorr_zero_pos {f : ℝ → ℝ} (hf : Continuous f) (hc : HasCompactSupport f)
    (hf0 : ∀ x, 0 ≤ f x) {x₀ : ℝ} (hx : f x₀ ≠ 0) : 0 < autocorr f 0 := by
  unfold autocorr
  simp only [sub_zero]
  exact (hf.mul hf).integral_pos_of_hasCompactSupport_nonneg_nonzero hc.mul_left
    (fun y => mul_nonneg (hf0 y) (hf0 y)) (x := x₀) (mul_ne_zero hx hx)

theorem halfWeight_prime_pos {p : ℕ} (hp : p.Prime) : 0 < halfWeight p := by
  unfold halfWeight
  rw [ArithmeticFunction.vonMangoldt_apply_prime hp]
  exact mul_pos (Real.log_pos (by exact_mod_cast hp.one_lt))
    (Real.rpow_pos_of_pos (by exact_mod_cast hp.pos) _)

/-- `stratumSN(log p) ≥ halfWeight(p)·‖f‖² > 0` for `f ≥ 0`, `f ≠ 0`. PROVED. -/
theorem stratumSN_pos (S : Finset ℕ) (N : ℕ) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S) (hpN : p ≤ N)
    {f : ℝ → ℝ} (hf : Continuous f) (hc : HasCompactSupport f) (hf0 : ∀ x, 0 ≤ f x)
    {x₀ : ℝ} (hx : f x₀ ≠ 0) :
    0 < stratumSN S N (autocorr f) (Real.log p) := by
  unfold stratumSN
  have hmem : (p, 1) ∈ cutIndex S N := by
    simp only [cutIndex, Finset.mem_filter, Finset.mem_product, Finset.mem_Icc, pow_one]
    exact ⟨⟨⟨hpS, hp⟩, le_refl 1, hp.one_lt.le.trans hpN⟩, hpN⟩
  have hterm : ∀ c ∈ cutIndex S N, 0 ≤ halfWeight (c.1 ^ c.2) *
      autocorr f (Real.log p - Real.log ((c.1 ^ c.2 : ℕ) : ℝ)) :=
    fun c _ => mul_nonneg (halfWeight_nonneg _) (autocorr_nonneg hf0 _)
  refine lt_of_lt_of_le ?_ (Finset.single_le_sum hterm hmem)
  show 0 < halfWeight (p ^ 1) * autocorr f (Real.log p - Real.log ((p ^ 1 : ℕ) : ℝ))
  rw [pow_one, sub_self]
  exact mul_pos (halfWeight_prime_pos hp) (autocorr_zero_pos hf hc hf0 hx)

/-- A smooth, compactly supported, nonnegative, nonzero test exists. PROVED. -/
theorem exists_smooth_test :
    ∃ f : ℝ → ℝ, ContDiff ℝ ∞ f ∧ HasCompactSupport f ∧ (∀ x, 0 ≤ f x) ∧ f 0 ≠ 0 := by
  let b : ContDiffBump (0 : ℝ) := ⟨1, 2, one_pos, one_lt_two⟩
  refine ⟨b, b.contDiff, b.hasCompactSupport, fun x => b.nonneg, ?_⟩
  rw [b.one_of_mem_closedBall (Metric.mem_closedBall_self zero_le_one)]
  norm_num

/-- The autocorrelation of a compactly supported function vanishes beyond some `R`. PROVED. -/
theorem autocorr_eventually_zero {f : ℝ → ℝ} (hc : HasCompactSupport f) :
    ∃ R : ℝ, ∀ v, R < v → autocorr f v = 0 := by
  obtain ⟨r, hr⟩ := hc.isCompact.isBounded.subset_closedBall 0
  refine ⟨2 * r, fun v hv => ?_⟩
  unfold autocorr
  have : (fun y => f y * f (y - v)) = fun _ => 0 := by
    funext y
    by_cases h1 : f y = 0
    · simp [h1]
    by_cases h2 : f (y - v) = 0
    · simp [h2]
    exfalso
    have m1 := hr (subset_tsupport f h1)
    have m2 := hr (subset_tsupport f h2)
    rw [Metric.mem_closedBall, Real.dist_eq, sub_zero] at m1 m2
    obtain ⟨a1, a2⟩ := abs_le.mp m1
    obtain ⟨b1, b2⟩ := abs_le.mp m2
    linarith
  rw [this, integral_zero]

section Datum2

variable {H Gu Gp : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
  [NormedAddCommGroup Gu] [InnerProductSpace ℂ Gu]
  [NormedAddCommGroup Gp] [InnerProductSpace ℂ Gp]

/-- **No covariant prime boundary.** For every datum with finite-cut BT3 whose truncation contains
a prime `p ≤ N`, some smooth compactly supported test admits no covariant prime boundary along its
orbit. PROVED, unconditional. -/
theorem no_covariant_prime_boundary (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S) (hpN : p ≤ N) :
    ∃ f : ℝ → ℝ, ContDiff ℝ ∞ f ∧ HasCompactSupport f ∧
      ∀ (U : ℝ → Gp →ₗᵢ[ℂ] Gp) (χ₀ χ₁ : ℝ → ℝ), (∀ W : ℝ, 0 ≤ W → χ₀ W * χ₁ W ≠ 0) →
        ¬ OrbitCovariantPrime D (autocorr f) U χ₀ χ₁ := by
  obtain ⟨f, hf, hc, hf0, hf1⟩ := exists_smooth_test
  obtain ⟨R, hR⟩ := autocorr_eventually_zero hc
  refine ⟨f, hf, hc, fun U χ₀ χ₁ hχ hcov => ?_⟩
  have hz := covariant_prime_forces_stratum_zero D S N hP f hf hc R hR U χ₀ χ₁ hcov hχ
    (Real.log p) (Real.log_natCast_nonneg p)
  have hpos := stratumSN_pos S N hp hpS hpN hf.continuous hc hf0 hf1
  linarith

end Datum2

/-- **Local-trace model.** If orbit vectors may jump only at a finite fixed cut set `C` and the
dilation is a translation, the jump set `J` of `ψ_0` obeys `J + W ⊆ C` for all `W ≥ 0`; then
`J = ∅`. PROVED. -/
theorem translate_jumpset_empty (C : Finset ℝ) (J : Set ℝ)
    (hJ : ∀ W : ℝ, 0 ≤ W → ∀ y ∈ J, y + W ∈ C) : J = ∅ := by
  by_contra hne
  obtain ⟨y, hy⟩ := Set.nonempty_iff_ne_empty.mpr hne
  have hinf : (C : Set ℝ).Infinite :=
    Set.infinite_of_injective_forall_mem (f := fun n : ℕ => y + n)
      (fun a b h => by simpa using h)
      (fun n => Finset.mem_coe.mpr (hJ n (Nat.cast_nonneg n) y hy))
  exact hinf C.finite_toSet

end GppBTSigmaBoundaryPW

#print axioms GppBTSigmaBoundaryPW.boundary_carriage_automatic
#print axioms GppBTSigmaBoundaryPW.punit_unitBlock_covariant
#print axioms GppBTSigmaBoundaryPW.primeBlock_covariant
#print axioms GppBTSigmaBoundaryPW.covariant_prime_forces_stratum_zero
#print axioms GppBTSigmaBoundaryPW.autocorr_nonneg
#print axioms GppBTSigmaBoundaryPW.autocorr_zero_pos
#print axioms GppBTSigmaBoundaryPW.halfWeight_prime_pos
#print axioms GppBTSigmaBoundaryPW.stratumSN_pos
#print axioms GppBTSigmaBoundaryPW.exists_smooth_test
#print axioms GppBTSigmaBoundaryPW.autocorr_eventually_zero
#print axioms GppBTSigmaBoundaryPW.no_covariant_prime_boundary
#print axioms GppBTSigmaBoundaryPW.translate_jumpset_empty
