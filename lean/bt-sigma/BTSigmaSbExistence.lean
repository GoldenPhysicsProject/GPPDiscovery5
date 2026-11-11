import GppVerify.RiemannHypothesis.BTSigmaNonCov

/-!
# BT-Σ: S_b existence, narrowed form — a datum with (P-COV), BT3, (P-UNIT), `Γ_u` graph-bounded and `H`-unbounded

Scratch module (Heart of the Cosmos / Grok, GPPDiscovery5 seat) for Turbo's question on
GPP-bridge PR #11, comments 6807739509 and 6808051667: does an infinite-dimensional datum exist
with (P-COV), BT3 and (P-UNIT), with `Γ_u` graph-bounded and `H`-unbounded? In boundary language:
construct a boundary-carried `P_W` with the forced shape
`greenForm ψ_W = e^{W/2}(−I c) + I·stratumSN(W)`, correction compactly supported in `W`, or name
the step that fails.

**Notions used (all from the tree, none invented).** `BTSigmaDatum` (BTSigma: `T*` as a
`LinearPMap`, `Γ₀u, Γ₁u, Γ₀p, Γ₁p`, Green's identity summed over the two blocks, joint
surjectivity, `testVec`); (P-COV) = `PCov`, BT3 = `BTSigmaTraceSN` (finite cut), (P-UNIT) = `PUnit`
(BTSigmaLocalEuler); "graph-bounded" = `‖Γ_iu x‖ ≤ K(‖x‖ + ‖T* x‖)` as in
`graphBounded_unit_product_upper`; "`H`-unbounded" = the negation of the hypotheses of
`bounded_unit_obstruction`; `opForm`, `perturbation_forced` (BTSigmaCovPerturbation);
`PNonCov`, `OrbitIndepCovariantPrime` (BTSigmaNonCov).

**The construction.** `H = (ℂ ⊕ L²(ℝ)) ⊕ (L²(ℝ) ⊕ L²(ℝ))`, vectors `(z, α, a, b)`; `Gu = ℂ`,
`Gp = L²(ℝ)`; `S_W = id ⊕ τ_W ⊕ id ⊕ τ_W` with `(τ_W g)(x) = g(x − W)` (a unitary group);
`dom T* = {e^{t/2} α ∈ L¹}`, `T*(z, α, a, b) = (i ℓ(α), 0, 0, −a)` with the unbounded functional
`ℓ(α) = ∫ e^{t/2} α`; `Γ₀u = z`, `Γ₁u = i ℓ(α)`, `Γ₀p = a`, `Γ₁p = b`. Test vectors
`ψ_W = S_W (1, α_h, f, b_h)` with `ℓ(α_h) = c/2`, `h = f ⋆ f̃` and
`b_h = −(i/2) Σ_{cuts} halfWeight(p^k) τ_{−log p^k} f`. The unit growth `e^{W/2}` is the
character `ℓ(τ_W α) = e^{W/2} ℓ(α)` of an `H`-unbounded, graph-bounded functional; the prime bump is
the translation pairing `⟪f, τ_s f⟫ = (f ⋆ f̃)(s)`.

**Results (all PROVED, standard axioms).**
* `green_identity`, `surj_boundary`, `datum`: a `BTSigmaDatum` for every truncation `(S, N)`.
* `datum_pcov` (P-COV), `datum_bt3` (finite-cut BT3, every smooth compactly supported `f`),
  `datum_punit` (P-UNIT).
* `datum_graphBounded`: `‖Γ_iu x‖ ≤ 1·(‖x‖ + ‖T* x‖)`. `Γ1u_unbounded`,
  `datum_unit_not_bounded`: `Γ_u` is `H`-unbounded. `H_not_finiteDimensional`.
* `datum_orbitIndep`: the prime boundary is covariant for independent intertwiners `U_W = id`,
  `V_W = τ_W` (the evader setting); `datum_pnoncov`: (P-NONCOV) on every nonnegative nonzero test.
* `datum_green_orbit`: the forced shape `greenForm ψ_W = e^{W/2}(−I c) + I·stratumSN(W)`.
* `datum_decomp`, `opForm_Aop`, `opForm_Pop`, `perturbation_constructed`: the perturbation
  `P_W x = (0, 0, 0, −τ_{−W} a)`, bounded (`Pop_norm_le`, norm ≤ 1) and acting on prime-block
  coordinates only, realises `T* ψ_W = S_W(e^{W/2} A ψ_0 + P_W ψ_0)` with
  `form_A(ψ_0) = −I c`, `form_{P_W}(ψ_0) = primeBlock(W) = I·stratumSN(W)`, zero past `R + log N`.
  `datum_boundary_carriage`: it meets every hypothesis of `boundary_carriage_automatic`.
* `sb_exists_narrowed`: the bundled answer.
* Scope (also PROVED): `Tst_eq_zero_on_ker`, `boundary_kernel_not_dense`: `T = T*|ker Γ` is `0`
  and `ker Γ` is not dense, so `T*` is **not** the adjoint of a densely defined symmetric operator.

**Verdict: CONSTRUCTED at the level of the typed clauses.** The narrowed question has a positive
answer: (P-COV), BT3, (P-UNIT), `Γ_u` graph-bounded and `H`-unbounded, `H` infinite-dimensional,
the forced-shape `P_W` boundary-carried and compactly supported in `W`, all hold at once. So no kill
of S_b can come from these clauses; the growth `e^{W/2}` vs unitarity is absorbed exactly by the
`H`-unbounded `Γ₁u`. What this datum does **not** meet, and where any kill must now live:
1. genuine BT2 — `T*` the adjoint of a densely defined closed symmetric `T` (dense `ker Γ`);
   untyped in `BTSigmaDatum`, fails here (`boundary_kernel_not_dense`). OPEN.
2. (P-LOG) provenance — the weights `halfWeight(p^k)` and cuts `log p^k` are inserted into the test
   amplitudes `b_h` (amplitude insertion, cf. `amplitude_insertion_reproduces_weights`). OPEN.
3. (P-θ) reflection `J` — not addressed. OPEN.
4. (P-LIM) — paper step, not formalised: for `f ≥ 0`, `‖b_h‖² ≥ ¼ Σ halfWeight(p^k)² ‖f‖²`, which
   diverges as `(S, N)` exhaust (`Σ log²p / p = ∞`), so this family fails (P-LIM).
S_b existence for a genuine datum, (P-LOG), (P-θ), (P-LIM), (BT-Σ), ARP and RH remain **OPEN**.

No `sorry`, no `admit`, no custom axiom; no external theorem is cited.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo.
-/

namespace GppBTSigmaSbExistence

open GppBTSigmaLocalEuler GppBTSigma GppArithmeticReflection Complex MeasureTheory
open GppBTSigmaCovObstruction GppBTSigmaBoundaryPW GppBTSigmaNonCov GppBTSigmaCovPerturbation
open GppBTSigmaCovHalfType GppBTSigmaCovGrowth GppBTSigmaC1
open scoped ContDiff

noncomputable section

/-- `L²(ℝ; ℂ)`. -/
abbrev L := Lp ℂ 2 (volume : Measure ℝ)

/-- Translation `(τ_W g)(x) = g(x − W)` on `L²(ℝ)`. -/
def shift (W : ℝ) : L →ₗᵢ[ℂ] L :=
  Lp.compMeasurePreservingₗᵢ ℂ (fun x : ℝ => x - W) (measurePreserving_sub_right volume W)

lemma coeFn_shift (W : ℝ) (g : L) :
    (shift W g : ℝ → ℂ) =ᵐ[volume] fun x => g (x - W) :=
  Lp.coeFn_compMeasurePreserving g (measurePreserving_sub_right volume W)

theorem shift_add (a b : ℝ) (g : L) : shift (a + b) g = shift a (shift b g) := by
  apply Lp.ext
  have h3 := (measurePreserving_sub_right volume a).quasiMeasurePreserving.ae_eq_comp
    (coeFn_shift b g)
  filter_upwards [coeFn_shift (a + b) g, coeFn_shift a (shift b g), h3] with x hx1 hx2 hx3
  rw [hx1, hx2]
  simp only [Function.comp] at hx3
  rw [hx3]
  congr 1; ring

theorem shift_zero (g : L) : shift 0 g = g := by
  apply Lp.ext
  filter_upwards [coeFn_shift 0 g] with x hx
  rw [hx, sub_zero]

/-- Translation as a linear isometric equivalence. -/
def shiftE (W : ℝ) : L ≃ₗᵢ[ℂ] L :=
  { toLinearEquiv := LinearEquiv.ofLinearMap (shift W).toLinearMap (shift (-W)).toLinearMap
      (LinearMap.ext fun g => by
        simp only [LinearMap.coe_comp, Function.comp_apply, LinearIsometry.coe_toLinearMap,
          LinearMap.id_coe, id_eq]
        rw [← shift_add, add_neg_cancel, shift_zero])
      (LinearMap.ext fun g => by
        simp only [LinearMap.coe_comp, Function.comp_apply, LinearIsometry.coe_toLinearMap,
          LinearMap.id_coe, id_eq]
        rw [← shift_add, neg_add_cancel, shift_zero])
    norm_map' := (shift W).norm_map }

@[simp] lemma shiftE_apply (W : ℝ) (g : L) : shiftE W g = shift W g := rfl

/-- The exponential weight `e^{t/2}`. -/
def wt (t : ℝ) : ℂ := (Real.exp (t / 2) : ℂ)

/-- `g` is integrable against `e^{t/2}` (the domain condition of the unit channel). -/
def WInt (g : L) : Prop := Integrable (fun t => wt t * g t)

/-- The unbounded unit functional `ℓ(g) = ∫ e^{t/2} g(t) dt`. -/
def ell (g : L) : ℂ := ∫ t, wt t * g t

lemma wt_shift (W t : ℝ) : wt t = wt W * wt (t - W) := by
  simp only [wt, ← Complex.ofReal_mul, ← Real.exp_add]
  congr 2; ring

theorem wint_shift {g : L} (hg : WInt g) (W : ℝ) : WInt (shift W g) := by
  have h := (hg.comp_sub_right W).const_mul (wt W)
  refine h.congr ?_
  filter_upwards [coeFn_shift W g] with x hx
  rw [hx, wt_shift W x, mul_assoc]

/-- `ℓ(τ_W g) = e^{W/2} ℓ(g)`: the unit functional is a character of the dilation. PROVED. -/
theorem ell_shift (g : L) (W : ℝ) : ell (shift W g) = wt W * ell g := by
  unfold ell
  have e1 : (∫ t, wt t * (shift W g) t) = ∫ t, wt W * (wt (t - W) * g (t - W)) := by
    refine integral_congr_ae ?_
    filter_upwards [coeFn_shift W g] with x hx
    rw [hx, wt_shift W x, mul_assoc]
  rw [e1, integral_const_mul]
  congr 1
  exact integral_sub_right_eq_self (fun t => wt t * g t) W

theorem wint_add {g₁ g₂ : L} (h₁ : WInt g₁) (h₂ : WInt g₂) : WInt (g₁ + g₂) := by
  refine (h₁.add h₂).congr ?_
  filter_upwards [Lp.coeFn_add g₁ g₂] with x hx
  simp only [hx, Pi.add_apply, mul_add]

theorem wint_smul (c : ℂ) {g : L} (h : WInt g) : WInt (c • g) := by
  refine (h.const_mul c).congr ?_
  filter_upwards [Lp.coeFn_smul c g] with x hx
  simp only [hx, Pi.smul_apply, smul_eq_mul]; ring

theorem wint_zero : WInt (0 : L) := by
  refine (integrable_zero ℝ ℂ volume).congr ?_
  filter_upwards [Lp.coeFn_zero ℂ 2 (volume : Measure ℝ)] with x hx
  rw [hx]; simp

theorem ell_add {g₁ g₂ : L} (h₁ : WInt g₁) (h₂ : WInt g₂) : ell (g₁ + g₂) = ell g₁ + ell g₂ := by
  unfold ell
  rw [← integral_add h₁ h₂]
  refine integral_congr_ae ?_
  filter_upwards [Lp.coeFn_add g₁ g₂] with x hx
  simp only [hx, Pi.add_apply, mul_add]

theorem ell_smul (c : ℂ) (g : L) : ell (c • g) = c * ell g := by
  unfold ell
  rw [← integral_const_mul]
  refine integral_congr_ae ?_
  filter_upwards [Lp.coeFn_smul c g] with x hx
  simp only [hx, Pi.smul_apply, smul_eq_mul]; ring

/-- The box `1_[0,1]`. -/
def box : L := indicatorConstLp 2 (measurableSet_Icc : MeasurableSet (Set.Icc (0:ℝ) 1))
  (by simp) (1 : ℂ)

lemma box_ae : (box : ℝ → ℂ) =ᵐ[volume] (Set.Icc (0:ℝ) 1).indicator fun _ => (1 : ℂ) :=
  indicatorConstLp_coeFn

lemma wt_box_ae : (fun t => wt t * box t) =ᵐ[volume] (Set.Icc (0:ℝ) 1).indicator wt := by
  filter_upwards [box_ae] with x hx
  rw [hx]
  by_cases h : x ∈ Set.Icc (0:ℝ) 1
  · simp [Set.indicator_of_mem h]
  · simp [Set.indicator_of_notMem h]

lemma wt_continuous : Continuous wt := by
  unfold wt; fun_prop

theorem wint_box : WInt box :=
  ((wt_continuous.integrableOn_Icc).integrable_indicator measurableSet_Icc).congr wt_box_ae.symm

/-- `κ = ℓ(1_[0,1])` is a real number `≥ 1`, hence nonzero. PROVED. -/
theorem ell_box : ∃ r : ℝ, 1 ≤ r ∧ ell box = r := by
  refine ⟨∫ t in Set.Icc (0:ℝ) 1, Real.exp (t / 2), ?_, ?_⟩
  · have h := setIntegral_ge_of_const_le_real (f := fun t => Real.exp (t / 2)) (c := 1)
      (measurableSet_Icc : MeasurableSet (Set.Icc (0:ℝ) 1))
      (by simp : (volume : Measure ℝ) (Set.Icc (0:ℝ) 1) ≠ ⊤)
      (fun t ht => Real.one_le_exp (by linarith [ht.1]))
      (by fun_prop : Continuous fun t : ℝ => Real.exp (t / 2)).integrableOn_Icc
    simpa using h
  · unfold ell
    rw [integral_congr_ae wt_box_ae, integral_indicator measurableSet_Icc]
    unfold wt
    exact integral_ofReal

theorem ell_box_ne : ell box ≠ 0 := by
  obtain ⟨r, hr, e⟩ := ell_box
  rw [e]; exact_mod_cast (by linarith : r ≠ 0)


/-! ## The carrier `H = (ℂ ⊕ L²) ⊕ (L² ⊕ L²)` and the dilation group -/

/-- Unit carrier `ℂ ⊕ L²(ℝ)`. -/
abbrev Hu := WithLp 2 (ℂ × L)
/-- Prime carrier `L²(ℝ) ⊕ L²(ℝ)`. -/
abbrev Hp := WithLp 2 (L × L)
/-- Total carrier. -/
abbrev H := WithLp 2 (Hu × Hp)

/-- Build a vector of `H` from its four components. -/
def mk (z : ℂ) (α a b : L) : H := WithLp.toLp 2 (WithLp.toLp 2 (z, α), WithLp.toLp 2 (a, b))

@[simp] lemma mk_z (z : ℂ) (α a b : L) : (mk z α a b).fst.fst = z := rfl
@[simp] lemma mk_α (z : ℂ) (α a b : L) : (mk z α a b).fst.snd = α := rfl
@[simp] lemma mk_a (z : ℂ) (α a b : L) : (mk z α a b).snd.fst = a := rfl
@[simp] lemma mk_b (z : ℂ) (α a b : L) : (mk z α a b).snd.snd = b := rfl

theorem H_ext {x y : H} (h1 : x.fst.fst = y.fst.fst) (h2 : x.fst.snd = y.fst.snd)
    (h3 : x.snd.fst = y.snd.fst) (h4 : x.snd.snd = y.snd.snd) : x = y := by
  rw [WithLp.ext_iff]
  refine Prod.ext ?_ ?_
  · show x.fst = y.fst
    rw [WithLp.ext_iff]; exact Prod.ext h1 h2
  · show x.snd = y.snd
    rw [WithLp.ext_iff]; exact Prod.ext h3 h4

/-- The dilation group `S_W = id_ℂ ⊕ τ_W ⊕ id_{L²} ⊕ τ_W`, a unitary group on `H`. -/
def Sd (W : ℝ) : H ≃ₗᵢ[ℂ] H :=
  LinearIsometryEquiv.withLpProdCongr 2
    (LinearIsometryEquiv.withLpProdCongr 2 (LinearIsometryEquiv.refl ℂ ℂ) (shiftE W))
    (LinearIsometryEquiv.withLpProdCongr 2 (LinearIsometryEquiv.refl ℂ L) (shiftE W))

@[simp] lemma Sd_z (W : ℝ) (x : H) : (Sd W x).fst.fst = x.fst.fst := rfl
@[simp] lemma Sd_α (W : ℝ) (x : H) : (Sd W x).fst.snd = shift W x.fst.snd := rfl
@[simp] lemma Sd_a (W : ℝ) (x : H) : (Sd W x).snd.fst = x.snd.fst := rfl
@[simp] lemma Sd_b (W : ℝ) (x : H) : (Sd W x).snd.snd = shift W x.snd.snd := rfl

/-- Group law `S_{a+b} = S_a S_b`. PROVED. -/
theorem Sd_add (a b : ℝ) (x : H) : Sd (a + b) x = Sd a (Sd b x) :=
  H_ext rfl (shift_add a b _) rfl (shift_add a b _)

/-! ## The domain, `T*` and the boundary maps -/

/-- `dom T* = {x : e^{t/2}·α ∈ L¹}` (condition on the unit `L²` component only). -/
def dom : Submodule ℂ H where
  carrier := {x | WInt x.fst.snd}
  add_mem' := fun {x y} hx hy => by
    show WInt (x + y).fst.snd
    rw [WithLp.add_fst, WithLp.add_snd]; exact wint_add hx hy
  zero_mem' := by
    show WInt (0 : H).fst.snd
    rw [WithLp.zero_fst, WithLp.zero_snd]; exact wint_zero
  smul_mem' := fun c x hx => by
    show WInt (c • x).fst.snd
    rw [WithLp.smul_fst, WithLp.smul_snd]; exact wint_smul c hx

/-- The unit functional on the domain. -/
def ellD : dom →ₗ[ℂ] ℂ where
  toFun x := ell (x : H).fst.snd
  map_add' x y := by
    show ell ((x : H) + (y : H)).fst.snd = _
    rw [WithLp.add_fst, WithLp.add_snd]; exact ell_add x.2 y.2
  map_smul' c x := by
    show ell (c • (x : H)).fst.snd = _
    rw [WithLp.smul_fst, WithLp.smul_snd]; exact ell_smul c _

/-- `T*(z, α, a, b) = (i ℓ(α), 0, 0, −a)`. -/
def Tlin : dom →ₗ[ℂ] H where
  toFun x := mk (I * ellD x) 0 0 (-(x : H).snd.fst)
  map_add' x y := by
    refine H_ext ?_ ?_ ?_ ?_ <;>
      simp [WithLp.add_fst, WithLp.add_snd, mul_add, map_add]
    abel
  map_smul' c x := by
    refine H_ext ?_ ?_ ?_ ?_ <;>
      simp [WithLp.smul_fst, WithLp.smul_snd, map_smul, mul_left_comm]

/-- The maximal operator of the construction. -/
def Tst : H →ₗ.[ℂ] H := ⟨dom, Tlin⟩

/-- `Γ₀u x = z` (bounded). -/
def Γ0u : Tst.domain →ₗ[ℂ] ℂ :=
  (WithLp.fstₗ 2 ℂ ℂ L) ∘ₗ (WithLp.fstₗ 2 ℂ Hu Hp) ∘ₗ dom.subtype
/-- `Γ₁u x = i ℓ(α)` (graph-bounded, `H`-unbounded). -/
def Γ1u : Tst.domain →ₗ[ℂ] ℂ := I • ellD
/-- `Γ₀p x = a`. -/
def Γ0p : Tst.domain →ₗ[ℂ] L :=
  (WithLp.fstₗ 2 ℂ L L) ∘ₗ (WithLp.sndₗ 2 ℂ Hu Hp) ∘ₗ dom.subtype
/-- `Γ₁p x = b`. -/
def Γ1p : Tst.domain →ₗ[ℂ] L :=
  (WithLp.sndₗ 2 ℂ L L) ∘ₗ (WithLp.sndₗ 2 ℂ Hu Hp) ∘ₗ dom.subtype

@[simp] lemma Γ0u_apply (x : Tst.domain) : Γ0u x = (x : H).fst.fst := rfl
@[simp] lemma Γ1u_apply (x : Tst.domain) : Γ1u x = I * ell (x : H).fst.snd := rfl
@[simp] lemma Γ0p_apply (x : Tst.domain) : Γ0p x = (x : H).snd.fst := rfl
@[simp] lemma Γ1p_apply (x : Tst.domain) : Γ1p x = (x : H).snd.snd := rfl
lemma Tst_apply (x : Tst.domain) :
    Tst x = mk (I * ell (x : H).fst.snd) 0 0 (-(x : H).snd.fst) := rfl

/-- **Green's identity** for the construction, on all of `dom T* × dom T*`. PROVED. -/
theorem green_identity (x y : Tst.domain) :
    greenForm Tst x y = blockForm (Γ0u x) (Γ1u x) (Γ0u y) (Γ1u y) +
      blockForm (Γ0p x) (Γ1p x) (Γ0p y) (Γ1p y) := by
  unfold greenForm blockForm
  rw [Tst_apply, Tst_apply]
  simp only [WithLp.prod_inner_apply, WithLp.ofLp_fst, WithLp.ofLp_snd, mk_z, mk_α, mk_a, mk_b,
    inner_zero_left, inner_zero_right, inner_neg_left, inner_neg_right, Γ0u_apply, Γ1u_apply,
    Γ0p_apply, Γ1p_apply, RCLike.inner_apply']
  ring

/-- **Surjectivity** of the boundary map onto `(ℂ × L²) × (ℂ × L²)`. PROVED. -/
theorem surj_boundary :
    Function.Surjective fun x : Tst.domain => ((Γ0u x, Γ0p x), (Γ1u x, Γ1p x)) := by
  rintro ⟨⟨u0, p0⟩, ⟨u1, p1⟩⟩
  set c : ℂ := -I * u1 / ell box
  have hm : mk u0 (c • box) p0 p1 ∈ dom := wint_smul c wint_box
  refine ⟨⟨mk u0 (c • box) p0 p1, hm⟩, Prod.ext (Prod.ext rfl rfl) (Prod.ext ?_ rfl)⟩
  show I * ell (c • box) = u1
  rw [ell_smul]
  have := ell_box_ne
  simp only [c]
  field_simp
  rw [I_sq]; ring

/-! ## Test vectors -/

open scoped Classical in
/-- `L²` image of a real continuous compactly supported function (`0` otherwise). -/
def toL2 (f : ℝ → ℝ) : L :=
  if h : Continuous f ∧ HasCompactSupport f then
    (Continuous.memLp_of_hasCompactSupport (continuous_ofReal.comp h.1)
      (h.2.comp_left Complex.ofReal_zero) : MemLp (fun x => (f x : ℂ)) 2 volume).toLp _
  else 0

lemma toL2_ae {f : ℝ → ℝ} (hf : Continuous f) (hc : HasCompactSupport f) :
    (toL2 f : ℝ → ℂ) =ᵐ[volume] fun x => (f x : ℂ) := by
  classical
  rw [toL2, dif_pos ⟨hf, hc⟩]
  exact MemLp.coeFn_toLp _

/-- `⟪f, τ_s f⟫ = (f ⋆ f̃)(s)`: the translation pairing is the autocorrelation. PROVED. -/
theorem inner_toL2_shift {f : ℝ → ℝ} (hf : Continuous f) (hc : HasCompactSupport f) (s : ℝ) :
    inner ℂ (toL2 f) (shift s (toL2 f)) = (autocorr f s : ℂ) := by
  rw [MeasureTheory.L2.inner_def]
  have h3 := (measurePreserving_sub_right volume s).quasiMeasurePreserving.ae_eq_comp
    (toL2_ae hf hc)
  have e : (fun a => inner ℂ ((toL2 f) a) ((shift s (toL2 f)) a)) =ᵐ[volume]
      fun a => ((f a * f (a - s) : ℝ) : ℂ) := by
    filter_upwards [toL2_ae hf hc, coeFn_shift s (toL2 f), h3] with a h1 h2 h3
    rw [h1, h2]
    simp only [Function.comp] at h3
    rw [h3, RCLike.inner_apply', Complex.conj_ofReal]
    push_cast; ring
  rw [integral_congr_ae e, integral_complex_ofReal]
  rfl

/-- `h` is the autocorrelation of a smooth compactly supported real function. -/
def IsAC (h : ℝ → ℝ) : Prop :=
  ∃ f : ℝ → ℝ, ContDiff ℝ ∞ f ∧ HasCompactSupport f ∧ autocorr f = h

open scoped Classical in
/-- A chosen square root `f` of `h = f ⋆ f̃` (any one; only `f ⋆ f̃ = h` is used). -/
def pick (h : ℝ → ℝ) : ℝ → ℝ := if hh : IsAC h then Classical.choose hh else 0

lemma pick_spec {f : ℝ → ℝ} (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f) :
    ContDiff ℝ ∞ (pick (autocorr f)) ∧ HasCompactSupport (pick (autocorr f)) ∧
      autocorr (pick (autocorr f)) = autocorr f := by
  classical
  have hh : IsAC (autocorr f) := ⟨f, hf, hc, rfl⟩
  rw [pick, dif_pos hh]
  exact Classical.choose_spec hh

/-- The cut position `log p^k`. -/
def cutL (c : ℕ × ℕ) : ℝ := Real.log ((c.1 ^ c.2 : ℕ) : ℝ)

/-- Prime channel, fixed component: `a_h = f`. -/
def ah (h : ℝ → ℝ) : L := toL2 (pick h)

/-- Prime channel, moving component: `b_h = −(i/2) Σ_cuts halfWeight(p^k) τ_{−log p^k} f`. -/
def bh (S : Finset ℕ) (N : ℕ) (h : ℝ → ℝ) : L :=
  (-(I / 2)) • ∑ c ∈ cutIndex S N, ((halfWeight (c.1 ^ c.2) : ℝ) : ℂ) • shift (-cutL c) (ah h)

/-- The unit integral `c_h = ∫ h(v) e^{−v/2} dv`. -/
def cU (h : ℝ → ℝ) : ℝ := ∫ v : ℝ, h v * Real.exp (-v / 2)

/-- Unit channel, moving component: `α_h = (c_h / 2κ)·1_[0,1]`, so `ℓ(α_h) = c_h/2`. -/
def αh (h : ℝ → ℝ) : L := (((cU h / 2 : ℝ) : ℂ) / ell box) • box

lemma ell_αh (h : ℝ → ℝ) : ell (αh h) = ((cU h / 2 : ℝ) : ℂ) := by
  rw [αh, ell_smul]; field_simp [ell_box_ne]

/-- The `W = 0` test vector `ψ_0 = (1, α_h, a_h, b_h)`. -/
def v0 (S : Finset ℕ) (N : ℕ) (h : ℝ → ℝ) : H := mk 1 (αh h) (ah h) (bh S N h)

lemma v0_mem (S : Finset ℕ) (N : ℕ) (h : ℝ → ℝ) : v0 S N h ∈ dom :=
  wint_smul _ wint_box

lemma Sd_mem {x : H} (hx : x ∈ dom) (W : ℝ) : Sd W x ∈ dom := by
  show WInt (Sd W x).fst.snd
  rw [Sd_α]; exact wint_shift hx W

/-- Test vectors `ψ_W = S_W ψ_0`. -/
def testVec (S : Finset ℕ) (N : ℕ) (h : ℝ → ℝ) (W : ℝ) : Tst.domain :=
  ⟨Sd W (v0 S N h), Sd_mem (v0_mem S N h) W⟩

/-! ## The datum -/

/-- **The datum.** `H = (ℂ ⊕ L²) ⊕ (L² ⊕ L²)`, `Gu = ℂ`, `Gp = L²(ℝ)`. -/
def datum (S : Finset ℕ) (N : ℕ) : BTSigmaDatum H ℂ L where
  Tstar := Tst
  Γ₀u := Γ0u
  Γ₁u := Γ1u
  Γ₀p := Γ0p
  Γ₁p := Γ1p
  green := green_identity
  surj := surj_boundary
  testVec := testVec S N

/-- **(P-COV)** holds with the unitary group `Sd`. PROVED. -/
theorem datum_pcov (S : Finset ℕ) (N : ℕ) : PCov (datum S N) Sd := by
  refine ⟨Sd_add, fun h W b => ?_⟩
  show Sd (W + b) (v0 S N h) = Sd b (Sd W (v0 S N h))
  rw [add_comm, Sd_add]

lemma tv_coe (S : Finset ℕ) (N : ℕ) (h : ℝ → ℝ) (W : ℝ) :
    (((datum S N).testVec h W : Tst.domain) : H) = Sd W (v0 S N h) := rfl

/-- **(P-UNIT)** holds (for every datum `h`, in particular every `f ⋆ f̃`). PROVED. -/
theorem datum_punit (S : Finset ℕ) (N : ℕ) : PUnit (datum S N) := by
  intro f _ _ W _
  show blockForm (Γ0u (testVec S N (autocorr f) W)) (Γ1u (testVec S N (autocorr f) W))
    (Γ0u (testVec S N (autocorr f) W)) (Γ1u (testVec S N (autocorr f) W)) = _
  have e0 : Γ0u (testVec S N (autocorr f) W) = 1 := rfl
  have e1 : Γ1u (testVec S N (autocorr f) W) =
      I * ((Real.exp (W / 2) * (cU (autocorr f) / 2) : ℝ) : ℂ) := by
    show I * ell (shift W (αh (autocorr f))) = _
    rw [ell_shift, ell_αh, wt]; push_cast; ring
  rw [e0, e1, unitStratumTerm]
  unfold blockForm
  simp only [RCLike.inner_apply', map_mul, Complex.conj_I, Complex.conj_ofReal, map_one]
  unfold cU
  push_cast; ring

/-- The prime pairing along the orbit: `⟪a_h, τ_W b_h⟫ = −(i/2)·stratumSN(W)`. PROVED. -/
theorem prime_pairing (S : Finset ℕ) (N : ℕ) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f)
    (hc : HasCompactSupport f) (W : ℝ) :
    inner ℂ (ah (autocorr f)) (shift W (bh S N (autocorr f))) =
      -(I / 2) * (stratumSN S N (autocorr f) W : ℂ) := by
  obtain ⟨hg, hgc, hga⟩ := pick_spec hf hc
  unfold bh
  rw [LinearIsometry.map_smul, inner_smul_right, map_sum, inner_sum]
  congr 1
  unfold stratumSN
  simp only [Complex.ofReal_sum, Complex.ofReal_mul]
  refine Finset.sum_congr rfl fun c _ => ?_
  rw [LinearIsometry.map_smul, inner_smul_right, ← shift_add, ah,
    inner_toL2_shift hg.continuous hgc, hga, cutL, ← sub_eq_add_neg]

/-- **Finite-cut BT3** holds: `primeBlock = I·stratumSN(W)` for every smooth compactly supported
`f` and `W ≥ 0` (indeed every real `W`). PROVED. -/
theorem datum_bt3 (S : Finset ℕ) (N : ℕ) : BTSigmaTraceSN (datum S N) S N := by
  intro f hf hc W _
  show blockForm (Γ0p (testVec S N (autocorr f) W)) (Γ1p (testVec S N (autocorr f) W))
    (Γ0p (testVec S N (autocorr f) W)) (Γ1p (testVec S N (autocorr f) W)) = _
  have e0 : Γ0p (testVec S N (autocorr f) W) = ah (autocorr f) := rfl
  have e1 : Γ1p (testVec S N (autocorr f) W) = shift W (bh S N (autocorr f)) := rfl
  rw [e0, e1]
  unfold blockForm
  rw [← inner_conj_symm, prime_pairing S N f hf hc W]
  simp only [map_mul, map_neg, map_div₀, Complex.conj_I, Complex.conj_ofReal, map_ofNat]
  ring

/-! ## Γ_u: graph-bounded, `H`-unbounded; `H` infinite-dimensional -/

/-- `‖Γ₀u x‖ ≤ ‖x‖`. PROVED. -/
theorem Γ0u_le (x : Tst.domain) : ‖Γ0u x‖ ≤ ‖(x : H)‖ :=
  (WithLp.norm_fst_le (x := (x : H).fst)).trans (WithLp.norm_fst_le (x := (x : H)))

/-- `‖Γ₁u x‖ ≤ ‖T* x‖`. PROVED. -/
theorem Γ1u_le (x : Tst.domain) : ‖Γ1u x‖ ≤ ‖Tst x‖ := by
  have e : Γ1u x = (Tst x).fst.fst := rfl
  rw [e]; exact (WithLp.norm_fst_le (x := (Tst x).fst)).trans (WithLp.norm_fst_le (x := Tst x))

/-- **`Γ_u` is bounded on the graph norm** (constant `1`). PROVED. -/
theorem datum_graphBounded (S : Finset ℕ) (N : ℕ) (x : (datum S N).Tstar.domain) :
    ‖(datum S N).Γ₀u x‖ ≤ 1 * (‖(x : H)‖ + ‖(datum S N).Tstar x‖) ∧
      ‖(datum S N).Γ₁u x‖ ≤ 1 * (‖(x : H)‖ + ‖(datum S N).Tstar x‖) := by
  refine ⟨?_, ?_⟩
  · have := Γ0u_le x
    have := norm_nonneg (Tst x)
    show ‖Γ0u x‖ ≤ 1 * (‖(x : H)‖ + ‖Tst x‖); linarith
  · have := Γ1u_le x
    have := norm_nonneg (x : H)
    show ‖Γ1u x‖ ≤ 1 * (‖(x : H)‖ + ‖Tst x‖); linarith

/-- **`Γ₁u` is unbounded in the `H`-norm.** Along `x_W = S_W(0, 1_[0,1], 0, 0)`, `‖x_W‖` is
constant while `|Γ₁u x_W| = e^{W/2} κ`. PROVED. -/
theorem Γ1u_unbounded : ¬ ∃ C : ℝ, ∀ x : Tst.domain, ‖Γ1u x‖ ≤ C * ‖(x : H)‖ := by
  rintro ⟨C, hC⟩
  obtain ⟨r, hr, er⟩ := ell_box
  have hm : mk 0 box 0 0 ∈ dom := wint_box
  set n := ‖mk 0 box 0 0‖
  set W : ℝ := 2 * (max C 0 * n)
  have hb := hC ⟨Sd W (mk 0 box 0 0), Sd_mem hm W⟩
  have hn : ‖((⟨Sd W (mk 0 box 0 0), Sd_mem hm W⟩ : Tst.domain) : H)‖ = n :=
    LinearIsometryEquiv.norm_map _ _
  have hv : Γ1u ⟨Sd W (mk 0 box 0 0), Sd_mem hm W⟩ = I * ((Real.exp (W / 2) * r : ℝ) : ℂ) := by
    show I * ell (shift W box) = _
    rw [ell_shift, er, wt]; push_cast; ring
  rw [hn, hv, norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Real.norm_eq_abs,
    abs_of_pos (by positivity)] at hb
  have hexp : W / 2 + 1 ≤ Real.exp (W / 2) := Real.add_one_le_exp _
  have hn0 : 0 ≤ n := norm_nonneg _
  have h1 : C * n ≤ max C 0 * n := mul_le_mul_of_nonneg_right (le_max_left _ _) hn0
  have h2 : Real.exp (W / 2) ≤ Real.exp (W / 2) * r :=
    le_mul_of_one_le_right (Real.exp_pos _).le hr
  have h3 : W / 2 = max C 0 * n := by simp only [W]; ring
  linarith

/-- Hence `Γ_u` is not `H`-bounded in the sense of `bounded_unit_obstruction`. PROVED. -/
theorem datum_unit_not_bounded (S : Finset ℕ) (N : ℕ) :
    ¬ ∃ C₀ C₁ : ℝ, (∀ x : (datum S N).Tstar.domain, ‖(datum S N).Γ₀u x‖ ≤ C₀ * ‖(x : H)‖) ∧
      (∀ x : (datum S N).Tstar.domain, ‖(datum S N).Γ₁u x‖ ≤ C₁ * ‖(x : H)‖) :=
  fun ⟨_, C₁, _, h₁⟩ => Γ1u_unbounded ⟨C₁, h₁⟩

/-- **`H` is infinite-dimensional** (forced: on a finite-dimensional space `Γ₁u` would be
bounded). PROVED. -/
theorem H_not_finiteDimensional : ¬ FiniteDimensional ℂ H := by
  intro hfd
  let Lc := LinearMap.toContinuousLinearMap Γ1u
  exact Γ1u_unbounded ⟨‖Lc‖, fun x => Lc.le_opNorm x⟩

/-! ## Non-covariance, the forced shape, and the perturbation `P_W` -/

/-- The prime boundary is covariant for **independent** intertwiners `U_W = id`, `V_W = τ_W`,
`χ₀ = χ₁ = 1` (the setting of `BTSigmaNonCov` / `BTSigmaQuasiAnalytic`). PROVED. -/
theorem datum_orbitIndep (S : Finset ℕ) (N : ℕ) (h : ℝ → ℝ) :
    OrbitIndepCovariantPrime (datum S N) h (fun _ => LinearIsometry.id) shift
      (fun _ => 1) (fun _ => 1) := by
  intro W _
  refine ⟨?_, ?_⟩
  · show (Sd W (v0 S N h)).snd.fst = ((1 : ℝ) : ℂ) • LinearIsometry.id ((Sd 0 (v0 S N h)).snd.fst)
    simp
  · show (Sd W (v0 S N h)).snd.snd = ((1 : ℝ) : ℂ) • shift W ((Sd 0 (v0 S N h)).snd.snd)
    simp [shift_zero]

/-- **(P-NONCOV)** holds on every nonnegative nonzero smooth test (via `pnoncov_of_pcov_bt3`):
no single isometry with a real nonvanishing character carries the prime block. PROVED. -/
theorem datum_pnoncov (S : Finset ℕ) (N : ℕ) {p : ℕ} (hp : p.Prime) (hpS : p ∈ S) (hpN : p ≤ N)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f) (hf0 : ∀ x, 0 ≤ f x)
    {x₀ : ℝ} (hx : f x₀ ≠ 0) : PNonCov (datum S N) Sd S N f :=
  pnoncov_of_pcov_bt3 (datum S N) Sd S N (datum_pcov S N) (datum_bt3 S N) hp hpS hpN f hf hc hf0 hx

/-- **The forced shape** on the orbit: `greenForm ψ_W = e^{W/2}(−I c) + I·stratumSN(W)`. PROVED. -/
theorem datum_green_orbit (S : Finset ℕ) (N : ℕ) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f)
    (hc : HasCompactSupport f) (W : ℝ) (hW : 0 ≤ W) :
    greenForm (datum S N).Tstar ((datum S N).testVec (autocorr f) W)
        ((datum S N).testVec (autocorr f) W) =
      (Real.exp (W / 2) : ℂ) * (-(I * (cU (autocorr f) : ℂ))) +
        I * (stratumSN S N (autocorr f) W : ℂ) := by
  rw [green_testVec_of_premises _ S N (datum_bt3 S N) (datum_punit S N) f hf hc W hW,
    unitStratumTerm, cU]
  push_cast; ring

/-- The fixed operator `A x = (κ z, 0, 0, 0)`. -/
def Aop (κ : ℂ) : H →ₗ[ℂ] H where
  toFun x := mk (κ * x.fst.fst) 0 0 0
  map_add' x y := by
    refine H_ext ?_ ?_ ?_ ?_ <;> simp [WithLp.add_fst, WithLp.add_snd, mul_add]
  map_smul' c x := by
    refine H_ext ?_ ?_ ?_ ?_ <;> simp [WithLp.smul_fst, WithLp.smul_snd, mul_left_comm]

/-- **The perturbation** `P_W x = (0, 0, 0, −τ_{−W} a)`: bounded, acting only on the prime-block
coordinates. -/
def Pop (W : ℝ) : H →ₗ[ℂ] H where
  toFun x := mk 0 0 0 (-(shift (-W) x.snd.fst))
  map_add' x y := by
    refine H_ext ?_ ?_ ?_ ?_ <;> simp [WithLp.add_fst, WithLp.add_snd, map_add]
    abel
  map_smul' c x := by
    refine H_ext ?_ ?_ ?_ ?_ <;> simp [WithLp.smul_fst, WithLp.smul_snd, map_smul]

/-- `‖P_W x‖ ≤ ‖x‖`: the perturbation is a bounded operator of norm `≤ 1`. PROVED. -/
theorem Pop_norm_le (W : ℝ) (x : H) : ‖Pop W x‖ ≤ ‖x‖ := by
  have e : ‖Pop W x‖ = ‖x.snd.fst‖ := by
    show ‖mk 0 0 0 (-(shift (-W) x.snd.fst))‖ = _
    have h1 := WithLp.prod_norm_sq_eq_of_L2 (mk 0 0 0 (-(shift (-W) x.snd.fst)))
    have h2 := WithLp.prod_norm_sq_eq_of_L2 (mk 0 0 0 (-(shift (-W) x.snd.fst))).fst
    have h3 := WithLp.prod_norm_sq_eq_of_L2 (mk 0 0 0 (-(shift (-W) x.snd.fst))).snd
    simp only [mk_z, mk_α, mk_a, mk_b, norm_zero, norm_neg, LinearIsometry.norm_map] at h2 h3
    have h4 : ‖mk 0 0 0 (-(shift (-W) x.snd.fst))‖ ^ 2 = ‖x.snd.fst‖ ^ 2 := by
      rw [h1, h2, h3]; ring
    have := norm_nonneg (mk 0 0 0 (-(shift (-W) x.snd.fst)))
    have := norm_nonneg x.snd.fst
    nlinarith [sq_nonneg (‖mk 0 0 0 (-(shift (-W) x.snd.fst))‖ - ‖x.snd.fst‖),
      sq_nonneg (‖mk 0 0 0 (-(shift (-W) x.snd.fst))‖ + ‖x.snd.fst‖)]
  rw [e]
  exact (WithLp.norm_fst_le (x := x.snd)).trans (WithLp.norm_snd_le (x := x))

/-- **Orbit decomposition** `T* ψ_W = S_W(e^{W/2} A ψ_0 + P_W ψ_0)` for every `W`. PROVED. -/
theorem datum_decomp (S : Finset ℕ) (N : ℕ) (h : ℝ → ℝ) (W : ℝ) :
    (datum S N).Tstar ((datum S N).testVec h W) =
      Sd W ((Real.exp (W / 2) : ℂ) • Aop (I * ((cU h / 2 : ℝ) : ℂ))
          ((datum S N).testVec h 0 : H) + Pop W ((datum S N).testVec h 0 : H)) := by
  show Tst (testVec S N h W) = _
  rw [Tst_apply]
  refine H_ext ?_ ?_ ?_ ?_
  · show I * ell (shift W (αh h)) = (Real.exp (W / 2) : ℂ) * (I * ((cU h / 2 : ℝ) : ℂ) * 1) + 0
    rw [ell_shift, ell_αh, wt]; ring
  · show (0 : L) = shift W ((Real.exp (W / 2) : ℂ) • (0 : L) + 0)
    simp
  · show (0 : L) = (Real.exp (W / 2) : ℂ) • (0 : L) + 0
    simp
  · show -(ah h) = shift W ((Real.exp (W / 2) : ℂ) • (0 : L) + -(shift (-W) (ah h)))
    rw [smul_zero, zero_add, LinearIsometry.map_neg, ← shift_add, add_neg_cancel, shift_zero]

/-- The fixed part has form `−I c` on `ψ_0`. PROVED. -/
theorem opForm_Aop (S : Finset ℕ) (N : ℕ) (h : ℝ → ℝ) :
    opForm (Aop (I * ((cU h / 2 : ℝ) : ℂ))) ((datum S N).testVec h 0 : H) = -(I * (cU h : ℂ)) := by
  show inner ℂ (mk (I * ((cU h / 2 : ℝ) : ℂ) * 1) 0 0 0) (Sd 0 (v0 S N h)) -
    inner ℂ (Sd 0 (v0 S N h)) (mk (I * ((cU h / 2 : ℝ) : ℂ) * 1) 0 0 0) = _
  simp only [WithLp.prod_inner_apply, WithLp.ofLp_fst, WithLp.ofLp_snd, mk_z, mk_α, mk_a, mk_b,
    Sd_z, inner_zero_left, inner_zero_right, RCLike.inner_apply', map_mul, Complex.conj_I,
    Complex.conj_ofReal, map_one, v0]
  push_cast; ring

lemma inner_shift_neg (W : ℝ) (a b : L) :
    inner ℂ (shift (-W) a) b = inner ℂ a (shift W b) := by
  rw [← (shift W).inner_map_map, ← shift_add, add_neg_cancel, shift_zero]

lemma inner_shift_neg' (W : ℝ) (a b : L) :
    inner ℂ a (shift (-W) b) = inner ℂ (shift W a) b := by
  rw [← (shift W).inner_map_map, ← shift_add, add_neg_cancel, shift_zero]

/-- **The perturbation form is the prime block**: `opForm(P_W)(ψ_0) = primeBlock(W)`,
for every `W`. PROVED. -/
theorem opForm_Pop (S : Finset ℕ) (N : ℕ) (h : ℝ → ℝ) (W : ℝ) :
    opForm (Pop W) ((datum S N).testVec h 0 : H) = primeBlock (datum S N) h W := by
  show inner ℂ (mk 0 0 0 (-(shift (-W) (ah h)))) (Sd 0 (v0 S N h)) -
    inner ℂ (Sd 0 (v0 S N h)) (mk 0 0 0 (-(shift (-W) (ah h)))) =
    inner ℂ (shift W (bh S N h)) (ah h) - inner ℂ (ah h) (shift W (bh S N h))
  simp only [WithLp.prod_inner_apply, WithLp.ofLp_fst, WithLp.ofLp_snd, mk_z, mk_α, mk_a, mk_b,
    Sd_z, Sd_α, Sd_a, Sd_b, inner_zero_left, inner_zero_right, inner_neg_left, inner_neg_right,
    v0, shift_zero, inner_shift_neg]
  rw [inner_shift_neg']
  ring

/-- **The perturbation, constructed.** For every smooth compactly supported `f`: the orbit
decomposition holds, the fixed part has form `−I c`, the perturbation form is `I·stratumSN(W)`
for `W ≥ 0`, and it vanishes past `R + log N` (compactly supported correction). PROVED. -/
theorem perturbation_constructed (S : Finset ℕ) (N : ℕ) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f)
    (hc : HasCompactSupport f) :
    (∀ W : ℝ, (datum S N).Tstar ((datum S N).testVec (autocorr f) W) =
      Sd W ((Real.exp (W / 2) : ℂ) • Aop (I * ((cU (autocorr f) / 2 : ℝ) : ℂ))
          ((datum S N).testVec (autocorr f) 0 : H) +
        Pop W ((datum S N).testVec (autocorr f) 0 : H))) ∧
    opForm (Aop (I * ((cU (autocorr f) / 2 : ℝ) : ℂ))) ((datum S N).testVec (autocorr f) 0 : H) =
      -(I * ((∫ v : ℝ, autocorr f v * Real.exp (-v / 2) : ℝ) : ℂ)) ∧
    (∀ W : ℝ, 0 ≤ W → opForm (Pop W) ((datum S N).testVec (autocorr f) 0 : H) =
      I * (stratumSN S N (autocorr f) W : ℂ)) ∧
    ∃ R : ℝ, ∀ W : ℝ, R + Real.log N < W →
      opForm (Pop W) ((datum S N).testVec (autocorr f) 0 : H) = 0 := by
  obtain ⟨R, hR⟩ := autocorr_eventually_zero hc
  refine ⟨datum_decomp S N _, opForm_Aop S N _, fun W hW => ?_, R, fun W hW => ?_⟩
  · rw [opForm_Pop]; exact datum_bt3 S N f hf hc W hW
  · rw [opForm_Pop]
    by_cases h0 : 0 ≤ W
    · rw [datum_bt3 S N f hf hc W h0]; exact perturbation_form_vanishes_late S N f R hR W hW
    · have hpb : primeBlock (datum S N) (autocorr f) W =
          I * (stratumSN S N (autocorr f) W : ℂ) := by
        show blockForm (ah (autocorr f)) (shift W (bh S N (autocorr f)))
          (ah (autocorr f)) (shift W (bh S N (autocorr f))) = _
        unfold blockForm
        rw [← inner_conj_symm, prime_pairing S N f hf hc W]
        simp only [map_mul, map_neg, map_div₀, Complex.conj_I, Complex.conj_ofReal, map_ofNat]
        ring
      rw [hpb]; exact perturbation_form_vanishes_late S N f R hR W hW

/-- Cross-check: the construction meets every hypothesis of `boundary_carriage_automatic`
(BTSigmaBoundaryPW), whose conclusion therefore holds for it. PROVED. -/
theorem datum_boundary_carriage (S : Finset ℕ) (N : ℕ) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f)
    (hc : HasCompactSupport f) :
    ∀ W : ℝ, 0 ≤ W →
      opForm (Pop W) ((datum S N).testVec (autocorr f) 0 : H) =
          primeBlock (datum S N) (autocorr f) W ∧
      (Real.exp (W / 2) : ℂ) *
          opForm (Aop (I * ((cU (autocorr f) / 2 : ℝ) : ℂ)))
            ((datum S N).testVec (autocorr f) 0 : H) =
        unitBlock (datum S N) (autocorr f) W := by
  obtain ⟨hdec, -, -, R, hlate⟩ := perturbation_constructed S N f hf hc
  obtain ⟨R', hR'⟩ := autocorr_eventually_zero hc
  exact boundary_carriage_automatic (datum S N) S N (datum_bt3 S N) (datum_punit S N) Sd
    (datum_pcov S N) f hf hc R' hR' _ Pop (fun W _ => hdec W) (R + Real.log N + 1)
    (fun W hW => hlate W (by linarith))

/-! ## Scope: what the construction does NOT give (the untyped BT2 clause) -/

/-- On the kernel of the boundary map `T*` vanishes: the "minimal operator" `T = T*|ker Γ` is
`0`. PROVED. -/
theorem Tst_eq_zero_on_ker (x : Tst.domain) (h1 : Γ1u x = 0) (h2 : Γ0p x = 0) : Tst x = 0 := by
  rw [Tst_apply]
  have e1 : I * ell (x : H).fst.snd = 0 := h1
  have e2 : (x : H).snd.fst = 0 := h2
  rw [e1, e2, neg_zero]
  refine H_ext ?_ ?_ ?_ ?_ <;> simp [WithLp.zero_fst, WithLp.zero_snd]

/-- **`ker Γ` is not dense in `H`**: every `x` with `Γ₀u x = 0` stays at distance `≥ 1` from
`(1, 0, 0, 0)`. So `T*` is NOT the adjoint of a densely defined symmetric operator; the
construction meets the typed `BTSigmaDatum` clauses only. PROVED. -/
theorem boundary_kernel_not_dense (x : Tst.domain) (h : Γ0u x = 0) :
    1 ≤ ‖mk 1 0 0 0 - (x : H)‖ := by
  have e : (mk 1 0 0 0 - (x : H)).fst.fst = 1 := by
    rw [WithLp.sub_fst, WithLp.sub_fst]
    have h' : (x : H).fst.fst = 0 := h
    rw [h', mk_z, sub_zero]
  have := (WithLp.norm_fst_le (x := (mk 1 0 0 0 - (x : H)).fst)).trans
    (WithLp.norm_fst_le (x := mk 1 0 0 0 - (x : H)))
  rw [e, norm_one] at this
  exact this

/-! ## The verdict -/

/-- **S_b existence, narrowed form: CONSTRUCTED.** For every truncation `(S, N)` there is an
infinite-dimensional datum with (P-COV), finite-cut BT3 and (P-UNIT), whose `Γ_u` is bounded on
the graph norm and unbounded in the `H`-norm. PROVED. -/
theorem sb_exists_narrowed (S : Finset ℕ) (N : ℕ) :
    ∃ (D : BTSigmaDatum H ℂ L) (Sb : ℝ → H ≃ₗᵢ[ℂ] H),
      PCov D Sb ∧ BTSigmaTraceSN D S N ∧ PUnit D ∧
      (∀ x : D.Tstar.domain, ‖D.Γ₀u x‖ ≤ 1 * (‖(x : H)‖ + ‖D.Tstar x‖) ∧
        ‖D.Γ₁u x‖ ≤ 1 * (‖(x : H)‖ + ‖D.Tstar x‖)) ∧
      (¬ ∃ C₀ C₁ : ℝ, (∀ x : D.Tstar.domain, ‖D.Γ₀u x‖ ≤ C₀ * ‖(x : H)‖) ∧
        (∀ x : D.Tstar.domain, ‖D.Γ₁u x‖ ≤ C₁ * ‖(x : H)‖)) ∧
      ¬ FiniteDimensional ℂ H :=
  ⟨datum S N, Sd, datum_pcov S N, datum_bt3 S N, datum_punit S N, datum_graphBounded S N,
    datum_unit_not_bounded S N, H_not_finiteDimensional⟩

end

end GppBTSigmaSbExistence

#print axioms GppBTSigmaSbExistence.shift_add
#print axioms GppBTSigmaSbExistence.shift_zero
#print axioms GppBTSigmaSbExistence.ell_shift
#print axioms GppBTSigmaSbExistence.wint_shift
#print axioms GppBTSigmaSbExistence.ell_box
#print axioms GppBTSigmaSbExistence.ell_box_ne
#print axioms GppBTSigmaSbExistence.green_identity
#print axioms GppBTSigmaSbExistence.surj_boundary
#print axioms GppBTSigmaSbExistence.Sd_add
#print axioms GppBTSigmaSbExistence.inner_toL2_shift
#print axioms GppBTSigmaSbExistence.datum_pcov
#print axioms GppBTSigmaSbExistence.datum_punit
#print axioms GppBTSigmaSbExistence.prime_pairing
#print axioms GppBTSigmaSbExistence.datum_bt3
#print axioms GppBTSigmaSbExistence.Γ0u_le
#print axioms GppBTSigmaSbExistence.Γ1u_le
#print axioms GppBTSigmaSbExistence.datum_graphBounded
#print axioms GppBTSigmaSbExistence.Γ1u_unbounded
#print axioms GppBTSigmaSbExistence.datum_unit_not_bounded
#print axioms GppBTSigmaSbExistence.H_not_finiteDimensional
#print axioms GppBTSigmaSbExistence.datum_orbitIndep
#print axioms GppBTSigmaSbExistence.datum_pnoncov
#print axioms GppBTSigmaSbExistence.datum_green_orbit
#print axioms GppBTSigmaSbExistence.Pop_norm_le
#print axioms GppBTSigmaSbExistence.datum_decomp
#print axioms GppBTSigmaSbExistence.opForm_Aop
#print axioms GppBTSigmaSbExistence.opForm_Pop
#print axioms GppBTSigmaSbExistence.perturbation_constructed
#print axioms GppBTSigmaSbExistence.datum_boundary_carriage
#print axioms GppBTSigmaSbExistence.Tst_eq_zero_on_ker
#print axioms GppBTSigmaSbExistence.boundary_kernel_not_dense
#print axioms GppBTSigmaSbExistence.sb_exists_narrowed
