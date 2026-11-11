import GppVerify.RiemannHypothesis.BTSigmaCovObstruction

/-!
# BT-Σ: genuine BT2 + (P-θ) — the centre kill (Lean) and the operator-level doubling (paper)

Scratch module (Heart of the Cosmos / Grok, GPPDiscovery5 seat) for Turbo's question on
GPP-bridge PR #11, comment 6808790481: a datum with dense `ker Γ` **and** a unitary `J` with
`J T J = −T`, plus (P-COV), finite-cut BT3, (P-UNIT), `Γ_u` graph-bounded and `H`-unbounded.

**(A) Tree (P-θ) clause — KILLED (PROVED).** `PLogTheta` asks for a unitary involution `J`
with `T* J = −J T*` on the domain **and** `testVec h (−W) = J (testVec h W)`. An anticommuting
unitary flips the Green form (`greenForm_anticomm`), so the `J`-fixed centre vector
`testVec h 0` has zero Green form (`greenForm_eq_zero_of_fixed`). With BT3 and (P-UNIT) at
`W = 0` this forces `stratumSN(h,0) = unitStratumTerm(h,0)` for every `h = f ⋆ f̃`
(`ptheta_centre_balance`). The concrete bump `φ` (support radius `1/4`) has zero prime stratum
at `W = 0` (`stratumSN_bump_zero`: every cut sits at `−log p^k ≤ −log 2 < −1/2`) but unit stratum
`(∫ φ e^{−y/2})(∫ φ e^{u/2}) > 0` (`integral_autocorr_mul_exp`, `bump_exp_pos`). Hence
`¬ (PLogTheta D J ∧ PUnit (D S N))` for every family, every `J`, every cut
(`ptheta_punit_kill`). Failing step: the reflection's fixed point `W = 0` cannot carry the
nonzero unit pole term; it needs neither density of `ker Γ` nor (P-COV), nor any property of `T`.

**(B) Operator-level clause only — CONSTRUCTED on paper (Lean port GAP).** Take the genuine
BT2 datum of `BTSigmaGenuineBT2` (`T ≥ 0`, half-line Laplacians) and double it:
`H' = H ⊕ H`, `T' = T ⊕ (−T)`, `J(x₁,x₂) = (x₂,x₁)`, so `J T' J = −T'` and `T'` is not
semibounded. Boundary maps `Γ₀' = (Γ₀x₁, Γ₀x₂)`, `Γ₁' = (Γ₁x₁, −Γ₁x₂)` blockwise:
Green's identity holds (the second copy's form is the negative), `Γ'` is onto, `ker Γ' =
dom T ⊕ dom T` is dense, each `Γ'_u` component is graph-bounded and `H`-unbounded. Test vectors
`(ψ_W, 0)` and `S'_b = S_b ⊕ S_b` give (P-COV), finite-cut BT3 and (P-UNIT) unchanged. This `J`
does **not** act `W ↦ −W` on the test family, and by (A) no `J` can; it is ζ-blind.

Lean results (7) PROVED, standard axioms. No `sorry`, no custom axiom, no cited theorem.
(P-LOG), ARP, CGC, RH and (BT-Σ) remain OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo.
-/

namespace GppBTSigmaGenuineBT2Theta

open GppBTSigmaLocalEuler GppBTSigma GppArithmeticReflection Complex MeasureTheory
open scoped InnerProductSpace ContDiff

section Green

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℂ E]

/-- An anticommuting unitary flips the sign of the Green form. PROVED. -/
theorem greenForm_anticomm (T : E →ₗ.[ℂ] E) (J : E ≃ₗᵢ[ℂ] E) (x y x' y' : T.domain)
    (hx : (x' : E) = J x) (hTx : T x' = -J (T x))
    (hy : (y' : E) = J y) (hTy : T y' = -J (T y)) :
    greenForm T x' y' = -greenForm T x y := by
  simp only [greenForm, hx, hTx, hy, hTy, inner_neg_left, inner_neg_right,
    LinearIsometryEquiv.inner_map_map]
  ring

/-- A `J`-fixed vector of the domain whose image is again in the domain with
`T (J x) = −J (T x)` has zero diagonal Green form. PROVED. -/
theorem greenForm_eq_zero_of_fixed (T : E →ₗ.[ℂ] E) (J : E ≃ₗᵢ[ℂ] E) (x : T.domain)
    (hanti : ∃ y : T.domain, (y : E) = J x ∧ T y = -J (T x)) (hfix : (x : E) = J x) :
    greenForm T x x = 0 := by
  obtain ⟨y, hy, hTy⟩ := hanti
  have hyx : y = x := Subtype.ext (hy.trans hfix.symm)
  subst hyx
  have e := greenForm_anticomm T J y y y y hy hTy hy hTy
  have : (2 : ℂ) * greenForm T y y = 0 := by linear_combination e
  simpa using this

end Green

section Datum

variable {H Gu Gp : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
  [NormedAddCommGroup Gu] [InnerProductSpace ℂ Gu]
  [NormedAddCommGroup Gp] [InnerProductSpace ℂ Gp]

/-- **Centre balance.** The (P-θ) clauses of `PLogTheta` (anticommuting `J`, test family
reflected `W ↦ −W`) together with (P-UNIT) force, at the fixed point `W = 0`, the truncated
prime stratum to equal the unit stratum. PROVED. -/
theorem ptheta_centre_balance (D : Finset ℕ → ℕ → BTSigmaDatum H Gu Gp) (J : H ≃ₗᵢ[ℂ] H)
    (hθ : PLogTheta D J) (S : Finset ℕ) (N : ℕ) (hU : PUnit (D S N))
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f) :
    stratumSN S N (autocorr f) 0 = unitStratumTerm (autocorr f) 0 := by
  obtain ⟨hBT3, -, hanti, hrefl⟩ := hθ
  have hfix := hrefl S N (autocorr f) 0
  rw [neg_zero] at hfix
  have hz := greenForm_eq_zero_of_fixed (D S N).Tstar J ((D S N).testVec (autocorr f) 0)
    (hanti S N _) hfix
  rw [green_testVec_of_premises (D S N) S N (hBT3 S N) hU f hf hc 0 le_rfl] at hz
  have h2 := (mul_eq_zero.mp hz).resolve_left Complex.I_ne_zero
  have h3 : stratumSN S N (autocorr f) 0 - unitStratumTerm (autocorr f) 0 = 0 := by
    exact_mod_cast h2
  linarith

end Datum

section Integral

/-- Factorisation of the unit-stratum integral of an autocorrelation:
`∫ (f ⋆ f̃)(v) e^{−v/2} dv = (∫ f e^{−y/2}) (∫ f e^{u/2})`. PROVED. -/
theorem integral_autocorr_mul_exp (f : ℝ → ℝ) (hf : Continuous f) (hc : HasCompactSupport f) :
    ∫ v : ℝ, autocorr f v * Real.exp (-v / 2) =
      (∫ y : ℝ, f y * Real.exp (-y / 2)) * ∫ u : ℝ, f u * Real.exp (u / 2) := by
  set C := ∫ u : ℝ, f u * Real.exp (u / 2)
  have h1 : ∀ v, autocorr f v * Real.exp (-v / 2) =
      ∫ y : ℝ, f y * f (y - v) * Real.exp (-v / 2) := by
    intro v; simp only [autocorr]; rw [integral_mul_const]
  simp_rw [h1]
  have hint : Integrable (Function.uncurry fun v y : ℝ => f y * f (y - v) * Real.exp (-v / 2))
      (volume.prod volume) := by
    apply Continuous.integrable_of_hasCompactSupport
    · exact (by fun_prop : Continuous fun p : ℝ × ℝ => f p.2 * f (p.2 - p.1) * Real.exp (-p.1 / 2))
    · set K := tsupport f
      have hK : IsCompact K := hc.isCompact
      set M := (fun q : ℝ × ℝ => q.1 - q.2) '' (K ×ˢ K)
      have hM : IsCompact M := (hK.prod hK).image (by fun_prop)
      refine HasCompactSupport.intro (hM.prod hK) ?_
      rintro ⟨v, y⟩ hp
      simp only [Function.uncurry_apply_pair]
      by_cases hy : y ∈ K
      · have hv : v ∉ M := fun h => hp ⟨h, hy⟩
        have : f (y - v) = 0 := by
          by_contra hne
          exact hv ⟨(y, y - v), ⟨hy, subset_tsupport f hne⟩, by simp⟩
        simp [this]
      · simp [image_eq_zero_of_notMem_tsupport hy]
  rw [integral_integral_swap hint]
  have h2 : ∀ y : ℝ, ∫ v : ℝ, f y * f (y - v) * Real.exp (-v / 2) =
      f y * Real.exp (-y / 2) * C := by
    intro y
    have e : ∫ v : ℝ, f (y - v) * Real.exp (-v / 2) = Real.exp (-y / 2) * C := by
      have := integral_sub_left_eq_self (fun u : ℝ => f u * Real.exp (-(y - u) / 2)) volume y
      simp only [sub_sub_cancel] at this
      rw [this, ← integral_const_mul]
      congr 1; ext u
      rw [show -(y - u) / 2 = -y / 2 + u / 2 by ring, Real.exp_add]; ring
    simp_rw [mul_assoc]
    rw [integral_const_mul, e]
  simp_rw [h2]
  rw [integral_mul_const]

end Integral

section Kill

/-- The concrete bump: radius `1/8` plateau, support the ball of radius `1/4`. -/
noncomputable def bump : ContDiffBump (0 : ℝ) := ⟨1 / 8, 1 / 4, by norm_num, by norm_num⟩

theorem bump_zero_outside (y : ℝ) (hy : (1 / 2 : ℝ) / 2 < |y|) : (bump : ℝ → ℝ) y = 0 := by
  have : y ∉ Function.support (bump : ℝ → ℝ) := by
    rw [ContDiffBump.support_eq]; simp only [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs]
    show ¬ |y| < 1 / 4; linarith
  simpa using this

/-- The small-support bump sees no prime cut at the centre `W = 0`. PROVED. -/
theorem stratumSN_bump_zero (S : Finset ℕ) (N : ℕ) :
    stratumSN S N (autocorr (bump : ℝ → ℝ)) 0 = 0 := by
  have hsupp := autocorr_hasSupportIn (f := (bump : ℝ → ℝ)) (L := 1 / 2) bump_zero_outside
  unfold stratumSN
  refine Finset.sum_eq_zero fun c hc => ?_
  simp only [cutIndex, Finset.mem_filter, Finset.mem_product, Finset.mem_Icc] at hc
  obtain ⟨⟨⟨_, hp⟩, hk1, _⟩, _⟩ := hc
  have h2 : (2 : ℝ) ≤ ((c.1 ^ c.2 : ℕ) : ℝ) := by
    have : 2 ≤ c.1 ^ c.2 := le_trans hp.two_le (Nat.le_self_pow (by omega) _)
    exact_mod_cast this
  have hl : Real.log 2 ≤ Real.log ((c.1 ^ c.2 : ℕ) : ℝ) := Real.log_le_log (by norm_num) h2
  have hl2 := Real.log_two_gt_d9
  rw [hsupp _ (by rw [zero_sub, abs_neg, abs_of_pos (by linarith)]; linarith), mul_zero]

theorem bump_exp_pos (a : ℝ) : 0 < ∫ y : ℝ, (bump : ℝ → ℝ) y * Real.exp (a * y) := by
  have hcont : Continuous fun y : ℝ => (bump : ℝ → ℝ) y * Real.exp (a * y) :=
    bump.continuous.mul (by fun_prop)
  refine hcont.integral_pos_of_hasCompactSupport_nonneg_nonzero (x := 0)
    (bump.hasCompactSupport.mul_right) (fun y => mul_nonneg bump.nonneg (Real.exp_pos _).le) ?_
  rw [bump.one_of_mem_closedBall (Metric.mem_closedBall_self (by simp [bump]))]
  simp

/-- **(P-θ) KILLED against (P-UNIT) + finite-cut BT3.** No family of genuine boundary triplets
carries the `PLogTheta` reflection (anticommuting unitary involution acting `W ↦ −W` on the test
family) together with (P-UNIT) at any cut: at the centre `W = 0` the small-support bump has
zero prime stratum but unit stratum `(∫ φ e^{−y/2})(∫ φ e^{u/2}) > 0`. PROVED. -/
theorem ptheta_punit_kill {H Gu Gp : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
    [NormedAddCommGroup Gu] [InnerProductSpace ℂ Gu]
    [NormedAddCommGroup Gp] [InnerProductSpace ℂ Gp]
    (D : Finset ℕ → ℕ → BTSigmaDatum H Gu Gp) (J : H ≃ₗᵢ[ℂ] H) (S : Finset ℕ) (N : ℕ) :
    ¬ (PLogTheta D J ∧ PUnit (D S N)) := by
  rintro ⟨hθ, hU⟩
  have e := ptheta_centre_balance D J hθ S N hU (bump : ℝ → ℝ) bump.contDiff bump.hasCompactSupport
  rw [stratumSN_bump_zero, unitStratumTerm,
    integral_autocorr_mul_exp _ bump.continuous bump.hasCompactSupport] at e
  have p1 := bump_exp_pos (-1 / 2)
  have p2 := bump_exp_pos (1 / 2)
  have q1 : (∫ y : ℝ, (bump : ℝ → ℝ) y * Real.exp (-y / 2)) =
      ∫ y : ℝ, (bump : ℝ → ℝ) y * Real.exp (-1 / 2 * y) := by
    congr 1; ext y; congr 2; ring
  have q2 : (∫ y : ℝ, (bump : ℝ → ℝ) y * Real.exp (y / 2)) =
      ∫ y : ℝ, (bump : ℝ → ℝ) y * Real.exp (1 / 2 * y) := by
    congr 1; ext y; congr 2; ring
  rw [q1, q2] at e
  have : 0 < Real.exp (0 / 2) * ((∫ y : ℝ, (bump : ℝ → ℝ) y * Real.exp (-1 / 2 * y)) *
      ∫ y : ℝ, (bump : ℝ → ℝ) y * Real.exp (1 / 2 * y)) :=
    mul_pos (Real.exp_pos _) (mul_pos p1 p2)
  linarith

end Kill

#print axioms greenForm_anticomm
#print axioms greenForm_eq_zero_of_fixed
#print axioms ptheta_centre_balance
#print axioms integral_autocorr_mul_exp
#print axioms stratumSN_bump_zero
#print axioms bump_exp_pos
#print axioms ptheta_punit_kill

end GppBTSigmaGenuineBT2Theta
