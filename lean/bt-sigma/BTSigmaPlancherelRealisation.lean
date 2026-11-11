/-
BTSigmaPlancherelRealisation.lean — GPPDiscovery5 (Heart of the Cosmos / Grok seat)
Task: Turbo, GPP-bridge PR #11 comment 6805592255, "Next exact question — the Plancherel
realisation": for an arbitrary C_c^∞ profile G, construct a, b ∈ L²(ℝ) with ⟪a, V_W b⟫ = G(W)
for all W, V_W the translation group, via the square-root (Paley–Wiener) factorisation.

Verdict: PROVED (all results standard axioms only, Mathlib-only imports, no sorry).

Exact characterisation of b̄a (Fourier side): with Mathlib's 𝓕 g(ξ) = ∫ 𝐞(−xξ) g(x) dx and
𝐞(t) = exp(2πit), a profile G is realisable as W ↦ ⟪a, V_W b⟫ iff G = 𝓕⁻¹h for some h ∈ L¹(ℝ),
and then conj(𝓕a)·(𝓕b) = h (the factorisation). For G continuous, integrable, 𝓕G integrable,
h = 𝓕G; choose 𝓕a = |h|^{1/2} (real, ≥ 0), 𝓕b = h/|h|^{1/2} (phase on b; carrying it on a instead
is the same up to conjugation). C_c^∞ ⊂ Schwartz ⊂ this class.

Contents:
* `mod W` : modulation (M_W g)(ξ) = 𝐞(ξW) g(ξ), a linear isometry of L²(ℝ); `mod_add` group law.
* `inner_mod_eq` : ⟪a, M_W b⟫ = 𝓕⁻¹(conj a · b)(W).
* `profile_is_fourier_L1` : necessity (conj a · b ∈ L¹ by Hölder).
* `realise_L1`, `realise_profile` : sufficiency, square-root factorisation.
* `V W := 𝓕⁻¹ ∘ M_W ∘ 𝓕` (Mathlib `Lp.fourierTransformₗᵢ`), `V_add`, `realise_profile_V`.
* `transl W` : (V_W b)(x) = b(x + W) (same definition as BTSigmaTranslationEvader).
  `transl_schwartz`, `fourier_transl_schwartz`, `fourier_transl` : 𝓕 ∘ transl W = M_W ∘ 𝓕 on L²
  (Schwartz density + `SchwartzMap.toLp_fourier_eq` + `fourierIntegral_comp_add_right`).
  `V_eq_transl` : V_W = transl W exactly.
* `realise_profile_transl`, `realise_schwartz_transl`, `realise_Cc_transl` : exact realisation
  ⟪a, transl W b⟫ = G(W) ∀W, for G continuous ∩ L¹ with 𝓕G ∈ L¹ / Schwartz / C_c^∞.
* `transl_profile_is_fourier_L1` : necessity through the translation group.
* `blockForm_profile` : for real σ (the stratumSN shape) with the integrability hypotheses,
  ∃ a b, −2i·Im⟪a, transl W b⟫ = i·σ(W) ∀W, i.e. blockForm(a, V_W b, a, V_W b) = i·stratumSN(W)
  via `blockForm_diag` (BTSigmaOneParam; cross-reference, not imported). U_W = id.

Not claimed: the quasi-analytic / Denjoy–Carleman follow-up (OPEN). No zeros used.
S_b existence, P-COV, P-UNIT, P-LOG, P-θ, P-LIM, BT-Σ, ARP and RH remain OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo.
-/
import Mathlib

open MeasureTheory Complex FourierTransform
open scoped Real ComplexInnerProductSpace ContDiff

noncomputable section

namespace BTSigmaPlancherel

/-- Gp = L²(ℝ; ℂ). -/
abbrev Gp := Lp ℂ 2 (volume : Measure ℝ)

/-- The character χ_W(ξ) = 𝐞(ξW) = exp(2πiξW). -/
def chi (W ξ : ℝ) : ℂ := (Real.fourierChar (ξ * W) : ℂ)

lemma norm_chi (W ξ : ℝ) : ‖chi W ξ‖ = 1 := by
  rw [chi, Real.fourierChar_apply]; exact Complex.norm_exp_ofReal_mul_I _

lemma continuous_chi (W : ℝ) : Continuous (chi W) := by
  unfold chi; fun_prop

lemma chi_add (W W' ξ : ℝ) : chi (W + W') ξ = chi W ξ * chi W' ξ := by
  simp [chi, mul_add, AddChar.map_add_eq_mul]

lemma memLp_chi_mul (W : ℝ) (f : Gp) : MemLp (fun ξ => chi W ξ * f ξ) 2 volume := by
  refine (Lp.memLp f).of_le ?_ (Filter.Eventually.of_forall fun ξ => ?_)
  · exact (continuous_chi W).aestronglyMeasurable.mul (Lp.aestronglyMeasurable f)
  · simp [norm_chi]

/-- Modulation (M_W f)(ξ) = χ_W(ξ) f(ξ), as a map on Gp. -/
def modFun (W : ℝ) (f : Gp) : Gp := (memLp_chi_mul W f).toLp _

lemma modFun_ae (W : ℝ) (f : Gp) : (modFun W f : ℝ → ℂ) =ᵐ[volume] fun ξ => chi W ξ * f ξ :=
  MemLp.coeFn_toLp _

/-- Modulation as a linear isometry. -/
def mod (W : ℝ) : Gp →ₗᵢ[ℂ] Gp where
  toFun := modFun W
  map_add' f g := by
    ext1
    filter_upwards [modFun_ae W (f + g), modFun_ae W f, modFun_ae W g, Lp.coeFn_add f g,
      Lp.coeFn_add (modFun W f) (modFun W g)] with ξ h1 h2 h3 h4 h5
    rw [h5, h1, h4, Pi.add_apply, Pi.add_apply, h2, h3]; ring
  map_smul' c f := by
    ext1
    filter_upwards [modFun_ae W (c • f), modFun_ae W f, Lp.coeFn_smul c f,
      Lp.coeFn_smul c (modFun W f)] with ξ h1 h2 h3 h4
    simp only [RingHom.id_apply]
    rw [h4, h1, h3]; simp only [Pi.smul_apply, smul_eq_mul, h2]; ring
  norm_map' f := by
    simp only [LinearMap.coe_mk, AddHom.coe_mk, modFun]
    rw [Lp.norm_toLp, Lp.norm_def]
    congr 1
    apply eLpNorm_congr_norm_ae
    filter_upwards with ξ
    simp [norm_chi]

lemma mod_ae (W : ℝ) (f : Gp) : (mod W f : ℝ → ℂ) =ᵐ[volume] fun ξ => chi W ξ * f ξ :=
  modFun_ae W f

/-- One-parameter group law M_{W+W'} = M_W ∘ M_{W'}. PROVED. -/
theorem mod_add (W W' : ℝ) (f : Gp) : mod (W + W') f = mod W (mod W' f) := by
  ext1
  filter_upwards [mod_ae (W + W') f, mod_ae W (mod W' f), mod_ae W' f] with ξ h1 h2 h3
  rw [h1, h2, h3, chi_add]; ring

/-- The pairing ⟪a, M_W b⟫ is the inverse Fourier transform of conj(a)·b. PROVED. -/
theorem inner_mod_eq (a b : Gp) (W : ℝ) :
    ⟪a, mod W b⟫ = 𝓕⁻ (fun ξ => (starRingEnd ℂ) (a ξ) * b ξ) W := by
  rw [L2.inner_def, Real.fourierInv_eq]
  apply integral_congr_ae
  filter_upwards [mod_ae W b] with ξ h
  rw [h]
  simp [chi, Circle.smul_def]
  ring_nf

/-- Necessity: every realised profile W ↦ ⟪a, M_W b⟫ is 𝓕⁻ of an L¹ function. PROVED. -/
theorem profile_is_fourier_L1 (a b : Gp) :
    ∃ h : ℝ → ℂ, Integrable h ∧ ∀ W, ⟪a, mod W b⟫ = 𝓕⁻ h W := by
  refine ⟨fun ξ => (starRingEnd ℂ) (a ξ) * b ξ, ?_, fun W => inner_mod_eq a b W⟩
  have ha : MemLp (fun ξ => (starRingEnd ℂ) (a ξ)) 2 volume := (Lp.memLp a).star
  exact ha.integrable_mul (Lp.memLp b)

/-! Sufficiency: the Paley–Wiener / square-root factorisation. -/

/-- |a| = |h|^{1/2}. -/
def sqrtAbs (h : ℝ → ℂ) (ξ : ℝ) : ℂ := (Real.sqrt ‖h ξ‖ : ℂ)
/-- b = h / |h|^{1/2} (the phase of h, times |h|^{1/2}; 0 where h = 0). -/
def phaseSqrt (h : ℝ → ℂ) (ξ : ℝ) : ℂ := h ξ / (Real.sqrt ‖h ξ‖ : ℂ)

lemma conj_sqrtAbs_mul_phaseSqrt (h : ℝ → ℂ) (ξ : ℝ) :
    (starRingEnd ℂ) (sqrtAbs h ξ) * phaseSqrt h ξ = h ξ := by
  unfold sqrtAbs phaseSqrt
  rw [Complex.conj_ofReal]
  by_cases h0 : h ξ = 0
  · simp [h0]
  · have : (Real.sqrt ‖h ξ‖ : ℂ) ≠ 0 := by
      have : 0 < Real.sqrt ‖h ξ‖ := Real.sqrt_pos.2 (norm_pos_iff.2 h0)
      exact_mod_cast this.ne'
    field_simp

lemma memLp_sqrtAbs {h : ℝ → ℂ} (hh : Integrable h) : MemLp (sqrtAbs h) 2 volume := by
  have hm : AEStronglyMeasurable (sqrtAbs h) volume := by
    unfold sqrtAbs
    exact (Complex.continuous_ofReal.comp_aestronglyMeasurable
      (hh.aestronglyMeasurable.norm.aemeasurable.sqrt.aestronglyMeasurable))
  rw [memLp_two_iff_integrable_sq_norm hm]
  refine hh.norm.congr (Filter.Eventually.of_forall fun ξ => ?_)
  simp [sqrtAbs, Real.sq_sqrt (norm_nonneg _)]

lemma memLp_phaseSqrt {h : ℝ → ℂ} (hh : Integrable h) : MemLp (phaseSqrt h) 2 volume := by
  refine (memLp_sqrtAbs hh).of_le ?_ (Filter.Eventually.of_forall fun ξ => ?_)
  · unfold phaseSqrt
    exact (hh.aestronglyMeasurable.aemeasurable.div
      (memLp_sqrtAbs hh).aestronglyMeasurable.aemeasurable).aestronglyMeasurable
  · unfold phaseSqrt sqrtAbs
    by_cases h0 : h ξ = 0
    · simp [h0]
    · have hs : 0 < Real.sqrt ‖h ξ‖ := Real.sqrt_pos.2 (norm_pos_iff.2 h0)
      rw [norm_div, Complex.norm_real, Real.norm_of_nonneg hs.le]
      rw [div_le_iff₀ hs, Real.mul_self_sqrt (norm_nonneg _)]

/-- Square-root factorisation: for every h ∈ L¹(ℝ) there are a, b ∈ L²(ℝ) with
⟪a, M_W b⟫ = 𝓕⁻ h (W) for all W. PROVED. -/
theorem realise_L1 {h : ℝ → ℂ} (hh : Integrable h) :
    ∃ a b : Gp, ∀ W, ⟪a, mod W b⟫ = 𝓕⁻ h W := by
  refine ⟨(memLp_sqrtAbs hh).toLp _, (memLp_phaseSqrt hh).toLp _, fun W => ?_⟩
  rw [inner_mod_eq]
  apply Real.fourierInv_congr_ae
  filter_upwards [MemLp.coeFn_toLp (memLp_sqrtAbs hh), MemLp.coeFn_toLp (memLp_phaseSqrt hh)]
    with ξ h1 h2
  rw [h1, h2, conj_sqrtAbs_mul_phaseSqrt]

/-- Exact-profile realisation: every continuous integrable G with integrable 𝓕 G (in particular
every C_c^∞ / Schwartz profile) is ⟪a, M_W b⟫ for some a, b ∈ L²(ℝ). PROVED. -/
theorem realise_profile {G : ℝ → ℂ} (hc : Continuous G) (hG : Integrable G)
    (hFG : Integrable (𝓕 G)) : ∃ a b : Gp, ∀ W, ⟪a, mod W b⟫ = G W := by
  obtain ⟨a, b, hab⟩ := realise_L1 hFG
  refine ⟨a, b, fun W => ?_⟩
  rw [hab, hc.fourierInv_fourier_eq hG hFG]

/-! Transfer through the L² Fourier transform (Mathlib `Lp.fourierTransformₗᵢ`). -/

/-- V_W := 𝓕⁻¹ ∘ M_W ∘ 𝓕 on Gp. -/
def V (W : ℝ) : Gp →ₗᵢ[ℂ] Gp :=
  ((Lp.fourierTransformₗᵢ ℝ ℂ).symm.toLinearIsometry.comp (mod W)).comp
    (Lp.fourierTransformₗᵢ ℝ ℂ).toLinearIsometry

/-- Group law for V. PROVED. -/
theorem V_add (W W' : ℝ) (f : Gp) : V (W + W') f = V W (V W' f) := by
  simp [V, mod_add]

/-- Realisation through V: same profile. PROVED. -/
theorem realise_profile_V {G : ℝ → ℂ} (hc : Continuous G) (hG : Integrable G)
    (hFG : Integrable (𝓕 G)) : ∃ a b : Gp, ∀ W, ⟪a, V W b⟫ = G W := by
  obtain ⟨a, b, hab⟩ := realise_profile hc hG hFG
  refine ⟨(Lp.fourierTransformₗᵢ ℝ ℂ).symm a, (Lp.fourierTransformₗᵢ ℝ ℂ).symm b, fun W => ?_⟩
  have : ⟪(Lp.fourierTransformₗᵢ ℝ ℂ).symm a, V W ((Lp.fourierTransformₗᵢ ℝ ℂ).symm b)⟫
      = ⟪a, mod W b⟫ := by
    simp only [V, LinearIsometry.coe_comp, Function.comp_apply,
      LinearIsometryEquiv.coe_toLinearIsometry, LinearIsometryEquiv.apply_symm_apply]
    exact (Lp.fourierTransformₗᵢ ℝ ℂ).symm.inner_map_map a (mod W b)
  rw [this, hab]


/-! Identification V_W = translation (Plancherel intertwining), via Schwartz density. -/
def transl (W : ℝ) : Gp →ₗᵢ[ℂ] Gp :=
  Lp.compMeasurePreservingₗᵢ ℂ (fun x : ℝ => x + W) (measurePreserving_add_right volume W)

lemma transl_schwartz (W : ℝ) (g : SchwartzMap ℝ ℂ) :
    transl W (g.toLp 2) = (SchwartzMap.compSubConstCLM ℂ (-W) g).toLp 2 := by
  show Lp.compMeasurePreserving (fun x : ℝ => x + W) (measurePreserving_add_right volume W)
    (g.toLp 2) = _
  ext1
  filter_upwards [Lp.coeFn_compMeasurePreserving (g.toLp 2) (measurePreserving_add_right volume W),
    SchwartzMap.coeFn_toLp (SchwartzMap.compSubConstCLM ℂ (-W) g) 2 volume,
    (measurePreserving_add_right volume W).quasiMeasurePreserving.ae
      (SchwartzMap.coeFn_toLp g 2 volume)] with x h1 h2 h3
  rw [h1, h2, Function.comp_apply, h3]
  simp [sub_neg_eq_add]

lemma fourier_transl_schwartz (W : ℝ) (g : SchwartzMap ℝ ℂ) (ξ : ℝ) :
    𝓕 (SchwartzMap.compSubConstCLM ℂ (-W) g) ξ = chi W ξ * 𝓕 g ξ := by
  rw [SchwartzMap.fourier_coe, SchwartzMap.fourier_coe]
  have : ((SchwartzMap.compSubConstCLM ℂ (-W) g : SchwartzMap ℝ ℂ) : ℝ → ℂ)
      = (g : ℝ → ℂ) ∘ fun v => v + W := by
    ext x; simp [sub_neg_eq_add]
  rw [this]
  have h := VectorFourier.fourierIntegral_comp_add_right Real.fourierChar volume
    (innerₗ ℝ) (g : ℝ → ℂ) W
  have h' := congrFun h ξ
  change VectorFourier.fourierIntegral _ _ _ _ ξ = chi W ξ * VectorFourier.fourierIntegral _ _ _ _ ξ
  rw [h', Circle.smul_def, smul_eq_mul]
  simp [chi, mul_comm]

/-- 𝓕 ∘ (translation by W) = M_W ∘ 𝓕 on L²(ℝ). PROVED (Schwartz density + Mathlib
`SchwartzMap.toLp_fourier_eq`, `fourierIntegral_comp_add_right`). -/
theorem fourier_transl (W : ℝ) (f : Gp) :
    𝓕 (transl W f) = mod W (𝓕 f) := by
  refine (SchwartzMap.denseRange_toLpCLM (E := ℝ) (F := ℂ) (p := 2) (μ := volume)
    ENNReal.ofNat_ne_top).induction_on f ?_ ?_
  · exact isClosed_eq ((Lp.fourierTransformₗᵢ ℝ ℂ).continuous.comp (transl W).continuous)
      ((mod W).continuous.comp (Lp.fourierTransformₗᵢ ℝ ℂ).continuous)
  · intro g
    simp only [SchwartzMap.toLpCLM_apply]
    rw [transl_schwartz, SchwartzMap.toLp_fourier_eq, SchwartzMap.toLp_fourier_eq]
    ext1
    filter_upwards [SchwartzMap.coeFn_toLp (𝓕 (SchwartzMap.compSubConstCLM ℂ (-W) g)) 2 volume,
      mod_ae W ((𝓕 g).toLp 2), SchwartzMap.coeFn_toLp (𝓕 g) 2 volume] with ξ h1 h2 h3
    rw [h1, h2, h3]
    exact fourier_transl_schwartz W g ξ

/-- V_W is exactly the translation group (V_W b)(x) = b(x + W). PROVED. -/
theorem V_eq_transl (W : ℝ) (f : Gp) : V W f = transl W f := by
  simp only [V, LinearIsometry.coe_comp, Function.comp_apply,
    LinearIsometryEquiv.coe_toLinearIsometry]
  have h := fourier_transl W f
  change (Lp.fourierTransformₗᵢ ℝ ℂ) (transl W f) = mod W ((Lp.fourierTransformₗᵢ ℝ ℂ) f) at h
  rw [← h, LinearIsometryEquiv.symm_apply_apply]

/-- Exact-profile realisation through the TRANSLATION group on L²(ℝ) (U_W = id, V_W = transl W):
for every continuous integrable G with integrable 𝓕 G (in particular C_c^∞ / Schwartz),
∃ a b ∈ L²(ℝ), ⟪a, V_W b⟫ = G(W) for all W. PROVED. -/
theorem realise_profile_transl {G : ℝ → ℂ} (hc : Continuous G) (hG : Integrable G)
    (hFG : Integrable (𝓕 G)) : ∃ a b : Gp, ∀ W, ⟪a, transl W b⟫ = G W := by
  obtain ⟨a, b, hab⟩ := realise_profile_V hc hG hFG
  exact ⟨a, b, fun W => by rw [← V_eq_transl, hab]⟩

/-- Necessity through the translation group: every ⟪a, V_W b⟫ is 𝓕⁻ of an L¹ function
(namely conj(𝓕a)·𝓕b). PROVED. -/
theorem transl_profile_is_fourier_L1 (a b : Gp) :
    ∃ h : ℝ → ℂ, Integrable h ∧ ∀ W, ⟪a, transl W b⟫ = 𝓕⁻ h W := by
  obtain ⟨h, hh, hW⟩ := profile_is_fourier_L1 (Lp.fourierTransformₗᵢ ℝ ℂ a)
    (Lp.fourierTransformₗᵢ ℝ ℂ b)
  refine ⟨h, hh, fun W => ?_⟩
  rw [← hW, ← V_eq_transl]
  simp only [V, LinearIsometry.coe_comp, Function.comp_apply,
    LinearIsometryEquiv.coe_toLinearIsometry]
  rw [← (Lp.fourierTransformₗᵢ ℝ ℂ).inner_map_map, LinearIsometryEquiv.apply_symm_apply]

/-- Schwartz (hence every C_c^∞) profile G : realised exactly through translations. PROVED. -/
theorem realise_schwartz_transl (G : SchwartzMap ℝ ℂ) :
    ∃ a b : Gp, ∀ W, ⟪a, transl W b⟫ = G W := by
  refine realise_profile_transl G.continuous G.integrable ?_
  have := (𝓕 G).integrable (μ := volume)
  rwa [SchwartzMap.fourier_coe] at this

/-- C_c^∞ profile (the stratumSN shape class): realised exactly through translations. PROVED. -/
theorem realise_Cc_transl {G : ℝ → ℂ} (h₁ : HasCompactSupport G) (h₂ : ContDiff ℝ ∞ G) :
    ∃ a b : Gp, ∀ W, ⟪a, transl W b⟫ = G W :=
  realise_schwartz_transl (h₁.toSchwartzMap h₂)

/-- The prime-block form profile: for a real profile σ (the stratumSN shape), with
G = −(i/2)σ, Im⟪a, V_W b⟫ = −σ(W)/2, i.e. −2i·Im⟪a, V_W b⟫ = i·σ(W)
(= blockForm(a, V_W b, a, V_W b) by `blockForm_diag`, BTSigmaOneParam). PROVED. -/
theorem blockForm_profile {σ : ℝ → ℝ} (hc : Continuous σ) (hσ : Integrable σ)
    (hFσ : Integrable (𝓕 (fun W => (-(I / 2)) * (σ W : ℂ)))) :
    ∃ a b : Gp, ∀ W, -2 * I * ((⟪a, transl W b⟫).im : ℂ) = I * (σ W : ℂ) := by
  have hc' : Continuous (fun W => (-(I / 2)) * (σ W : ℂ)) := by fun_prop
  have hG : Integrable (fun W => (-(I / 2)) * (σ W : ℂ)) := (hσ.ofReal).const_mul _
  obtain ⟨a, b, hab⟩ := realise_profile_transl hc' hG hFσ
  refine ⟨a, b, fun W => ?_⟩
  rw [hab]
  simp
  ring

end BTSigmaPlancherel

#print axioms BTSigmaPlancherel.mod_add
#print axioms BTSigmaPlancherel.inner_mod_eq
#print axioms BTSigmaPlancherel.profile_is_fourier_L1
#print axioms BTSigmaPlancherel.realise_L1
#print axioms BTSigmaPlancherel.realise_profile
#print axioms BTSigmaPlancherel.V_add
#print axioms BTSigmaPlancherel.realise_profile_V
#print axioms BTSigmaPlancherel.transl_schwartz
#print axioms BTSigmaPlancherel.fourier_transl_schwartz
#print axioms BTSigmaPlancherel.fourier_transl
#print axioms BTSigmaPlancherel.V_eq_transl
#print axioms BTSigmaPlancherel.realise_profile_transl
#print axioms BTSigmaPlancherel.transl_profile_is_fourier_L1
#print axioms BTSigmaPlancherel.realise_schwartz_transl
#print axioms BTSigmaPlancherel.realise_Cc_transl
#print axioms BTSigmaPlancherel.blockForm_profile
