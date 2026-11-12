/-
Copyright (c) 2026 Daniel Toupin. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Daniel Toupin (with AI assistance: Heart of the Cosmos / Grok)
-/
import Mathlib

/-!
# `WeilData` is inhabited, and `HPsiEquiv` is a theorem (not a premise)

Scratchpad version of GPPVerify `GppVerify/RiemannHypothesis/KreinHPsiWeilData.lean`
(branch `heart/hpsi`). Self-contained: the definitions `kform`, `HAllWindows`, `window`,
`fit_translates` (GPPVerify KreinBridge) and `triPsi … HPsiEquiv` (GPPVerify KreinHPsi, 65358d6)
are copied verbatim below so this file compiles on its own.

Labels:
* `tau_linearIndependent` (distinct triangle translates are linearly independent): PROVED.
* `exists_translate_data` (any kernel / cosh / sinh values are realised on translates by some
  bilinear `B` and linear `ℓ_c`, `ℓ_s` on `ℝ → ℝ`): PROVED. Kernel-blind.
* `weilData_nonempty`, `weilData`: PROVED / CONSTRUCTED from `Kψ` (pole, primes, Γ) and the
  cosh/sinh values. No zero of ζ is an input. The witness `B` is fixed on `F_ψ` and arbitrary
  on a complement; only its `F_ψ` part enters `HPsiWindows`.
* `hPsiWindows_iff`, `hPsiEquiv`: PROVED. The former named premise `HPsiEquiv` is discharged.
* `exists_weilData_iff`: PROVED (non-vacuous form).
* `HPsi` (= H_ψ): OPEN. ARP: OPEN. RH: OPEN.
-/

namespace GppKreinBridge

/-- The Toeplitz form of a kernel on a finite family of points. -/
def kform (K : ℝ → ℝ) {n : ℕ} (x : Fin n → ℝ) (c : Fin n → ℝ) : ℝ :=
  ∑ i, ∑ j, c i * c j * K (x i - x j)


variable {F : Type*} [AddCommGroup F] [Module ℝ F]

/-- Turbo's claim `T ≥ 2ℓ_s²` on `ker ℓ_c`, at every window, with
`T = Q − 2ℓ_c² + 2ℓ_s²`. -/
def HAllWindows (B : F →ₗ[ℝ] F →ₗ[ℝ] ℝ) (S : ℝ → Submodule ℝ F)
    (lc ls : ℝ → F →ₗ[ℝ] ℝ) : Prop :=
  ∀ L : ℝ, ∀ v ∈ S L, lc L v = 0 →
    2 * (ls L v) ^ 2 ≤ B v v - 2 * (lc L v) ^ 2 + 2 * (ls L v) ^ 2


/-- The concrete window: functions vanishing outside `[−L, L]`. -/
def window (L : ℝ) : Submodule ℝ (ℝ → ℝ) where
  carrier := {f | ∀ t, L < |t| → f t = 0}
  zero_mem' := fun _ _ => rfl
  add_mem' := fun hf hg t ht => by simp [hf t ht, hg t ht]
  smul_mem' := fun c f hf t ht => by simp [hf t ht]

/-- **Every finite family fits in one window**: translates of a function supported in
`[−a, a]` by `x₁, …, xₙ` all lie in the window `L = a + Σ |xᵢ|`. -/
theorem fit_translates (ψ : ℝ → ℝ) (a : ℝ) (hψ : ∀ t, a < |t| → ψ t = 0) (n : ℕ)
    (x : Fin n → ℝ) : ∃ L, ∀ i, (fun t => ψ (t - x i)) ∈ window L := by
  refine ⟨a + ∑ j, |x j|, fun i t ht => hψ _ ?_⟩
  have h1 : |x i| ≤ ∑ j, |x j| :=
    Finset.single_le_sum (f := fun j => |x j|) (fun j _ => abs_nonneg _) (Finset.mem_univ i)
  have h2 : |t| ≤ |t - x i| + |x i| := by
    calc |t| = |(t - x i) + x i| := by ring_nf
      _ ≤ |t - x i| + |x i| := abs_add_le _ _
  linarith

end GppKreinBridge

namespace GppKreinHPsi

open GppKreinBridge MeasureTheory

/-- The triangle test function `ψ(t) = max(0, 1 − |t|)`, supported in `[−1, 1]`. -/
noncomputable def triPsi (t : ℝ) : ℝ := max 0 (1 - |t|)

/-- Translate `τ_x ψ = ψ(· − x)`. -/
noncomputable def tau (x : ℝ) : ℝ → ℝ := fun t => triPsi (t - x)

/-- Fourier transform of the triangle on the real axis, `Ψ(t) = (sin(t/2)/(t/2))²`. -/
noncomputable def bigPsi (t : ℝ) : ℝ := if t = 0 then 1 else (Real.sin (t / 2) / (t / 2)) ^ 2

/-- `Ψ(i/2) = Ψ(−i/2) = 16 sinh²(1/4)`. -/
noncomputable def bigPsiI2 : ℝ := 16 * Real.sinh (1 / 4) ^ 2

/-- Autocorrelation `g = ψ ⋆ ψ` in closed form. -/
noncomputable def triG (y : ℝ) : ℝ :=
  if |y| ≤ 1 then 2 / 3 - y ^ 2 + |y| ^ 3 / 2
  else if |y| ≤ 2 then (2 - |y|) ^ 3 / 6 else 0

/-- Pole term `2 Ψ(i/2)² cosh(d/2)`. -/
noncomputable def poleTerm (d : ℝ) : ℝ := 2 * bigPsiI2 ^ 2 * Real.cosh (d / 2)

/-- Prime term `Σ_n Λ(n) n^{−1/2} [g(log n + d) + g(−log n + d)]` (finitely many nonzero terms,
since `g` vanishes off `[−2, 2]`). -/
noncomputable def primeTerm (d : ℝ) : ℝ :=
  ∑' n : ℕ, ArithmeticFunction.vonMangoldt n / Real.sqrt n *
    (triG (Real.log n + d) + triG (-Real.log n + d))

/-- Archimedean term `(1/π) ∫_0^∞ Ψ(t)² cos(t d) (Re ψ_Γ(1/4 + it/2) − log π) dt`. -/
noncomputable def archTerm (d : ℝ) : ℝ :=
  (1 / Real.pi) * ∫ t in Set.Ioi (0 : ℝ), bigPsi t ^ 2 * Real.cos (t * d) *
    ((logDeriv Complex.Gamma ((1 / 4 : ℂ) + (t / 2 : ℝ) * Complex.I)).re - Real.log Real.pi)

/-- The translate kernel `Kψ(d) = B(τ_0 ψ, τ_d ψ)` on the zero-free side. -/
noncomputable def Kpsi (d : ℝ) : ℝ := poleTerm d - primeTerm d + archTerm d

/-- **H_ψ, Toeplitz form** (OPEN): positivity of `kform Kψ` on the cosh hyperplane. -/
def HPsi : Prop :=
  ∀ (n : ℕ) (x c : Fin n → ℝ), ∑ i, c i * Real.cosh (x i / 2) = 0 → 0 ≤ kform Kpsi x c

/-- `F_ψ = span_ℝ {τ_x ψ}`. -/
noncomputable def Fpsi : Submodule ℝ (ℝ → ℝ) := Submodule.span ℝ (Set.range tau)

/-- The windowed spaces `S L = F_ψ ∩ window L`. -/
noncomputable def Spsi (L : ℝ) : Submodule ℝ (ℝ → ℝ) := Fpsi ⊓ window L

/-- The data `HAllWindows` consumes on `F_ψ`: a bilinear `B` whose translate kernel is `Kψ`,
and L-independent `ℓ_c, ℓ_s` with the cosh / sinh pairings on translates. -/
structure WeilData where
  B : (ℝ → ℝ) →ₗ[ℝ] (ℝ → ℝ) →ₗ[ℝ] ℝ
  lc : (ℝ → ℝ) →ₗ[ℝ] ℝ
  ls : (ℝ → ℝ) →ₗ[ℝ] ℝ
  hK : ∀ x y, B (tau x) (tau y) = Kpsi (x - y)
  hlc : ∀ x, lc (tau x) = bigPsiI2 * Real.cosh (x / 2)
  hls : ∀ x, ls (tau x) = bigPsiI2 * Real.sinh (x / 2)

/-- **H_ψ, `HAllWindows` form** (OPEN) for given Weil data. -/
def HPsiWindows (W : WeilData) : Prop :=
  HAllWindows W.B Spsi (fun _ => W.lc) (fun _ => W.ls)

/-- **Named premise** (not proved here): the two forms agree for every choice of Weil data. -/
def HPsiEquiv : Prop := ∀ W : WeilData, HPsiWindows W ↔ HPsi

end GppKreinHPsi

namespace GppKreinHPsi

open GppKreinBridge

/-! ## 1. The triangle translates are linearly independent (PROVED) -/

theorem triPsi_eq_zero_of_one_le {t : ℝ} (h : 1 ≤ |t|) : triPsi t = 0 := by
  unfold triPsi; exact max_eq_left (by linarith)

theorem triPsi_of_abs_le {t : ℝ} (h : |t| ≤ 1) : triPsi t = 1 - |t| := by
  unfold triPsi; exact max_eq_right (by linarith)

/-- Distinct translates of `ψ` are linearly independent in `ℝ → ℝ`: evaluate just left of the
right end `a + 1` of the right-most translate, where every other translate already vanishes. -/
theorem tau_linearIndependent : LinearIndependent ℝ tau := by
  classical
  rw [linearIndependent_iff']
  intro s
  induction s using Finset.induction_on_max with
  | empty => intro g _ i hi; simp at hi
  | insert a s has ih =>
    intro g hg
    obtain ⟨δ, hδ0, hδ1, hδs⟩ : ∃ δ : ℝ, 0 < δ ∧ δ ≤ 1 ∧ ∀ x ∈ s, x ≤ a - δ := by
      rcases s.eq_empty_or_nonempty with h | h
      · exact ⟨1, one_pos, le_rfl, by simp [h]⟩
      · refine ⟨min 1 (a - s.max' h), lt_min one_pos (sub_pos.2 (has _ (s.max'_mem h))),
          min_le_left _ _, fun x hx => ?_⟩
        have h1 := s.le_max' x hx
        have h2 := min_le_right 1 (a - s.max' h)
        linarith
    have has' : a ∉ s := fun h => lt_irrefl a (has a h)
    rw [Finset.sum_insert has'] at hg
    have hev := congrFun hg (a + 1 - δ)
    simp only [Pi.add_apply, Finset.sum_apply, Pi.smul_apply, smul_eq_mul,
      Pi.zero_apply] at hev
    have h0 : ∀ x ∈ s, tau x (a + 1 - δ) = 0 := by
      intro x hx
      apply triPsi_eq_zero_of_one_le
      have := hδs x hx
      rw [abs_of_nonneg (by linarith)]; linarith
    rw [Finset.sum_eq_zero (fun x hx => by rw [h0 x hx, mul_zero])] at hev
    have hta : tau a (a + 1 - δ) = δ := by
      show triPsi (a + 1 - δ - a) = δ
      rw [triPsi_of_abs_le] <;> rw [abs_of_nonneg (by linarith)] <;> linarith
    rw [hta, add_zero] at hev
    have hga : g a = 0 := (mul_eq_zero.1 hev).resolve_right hδ0.ne'
    rw [hga, zero_smul, zero_add] at hg
    intro i hi
    rcases Finset.mem_insert.1 hi with rfl | hi
    · exact hga
    · exact ih g hg i hi

/-! ## 2. Prescribing linear and bilinear maps on a linearly independent family (PROVED) -/

section Prescribe

variable {V W : Type*} [AddCommGroup V] [Module ℝ V] [AddCommGroup W] [Module ℝ W] {ι : Type*}

/-- A linear map can take any prescribed values on a linearly independent family. -/
theorem exists_linear_on_li {v : ι → V} (hv : LinearIndependent ℝ v) (f : ι → W) :
    ∃ ℓ : V →ₗ[ℝ] W, ∀ i, ℓ (v i) = f i := by
  obtain ⟨g, hg⟩ := LinearMap.exists_extend ((Module.Basis.span hv).constr ℝ f)
  refine ⟨g, fun i => ?_⟩
  have h := LinearMap.congr_fun hg (Module.Basis.span hv i)
  rw [LinearMap.comp_apply, Module.Basis.constr_basis, Submodule.subtype_apply,
    Module.Basis.coe_span_apply] at h
  exact h

/-- A bilinear form can take any prescribed values on pairs from a linearly independent
family. -/
theorem exists_bilinear_on_li {v : ι → V} (hv : LinearIndependent ℝ v) (F : ι → ι → ℝ) :
    ∃ B : V →ₗ[ℝ] V →ₗ[ℝ] ℝ, ∀ i j, B (v i) (v j) = F i j := by
  choose row hrow using fun i => exists_linear_on_li hv (F i)
  obtain ⟨B, hB⟩ := exists_linear_on_li hv row
  exact ⟨B, fun i j => by rw [hB, hrow]⟩

end Prescribe

/-! ## 3. `WeilData` is inhabited (PROVED; the witness is CONSTRUCTED from `Kψ` alone) -/

/-- Kernel-blind form of the non-vacuity statement: for **any** kernel `K` and any
`ℓ_c, ℓ_s` values on translates there is a bilinear `B` and functionals with those values.
So inhabiting `WeilData` costs nothing beyond `tau_linearIndependent`, and carries no
positivity information. -/
theorem exists_translate_data (K fc fs : ℝ → ℝ) :
    ∃ (B : (ℝ → ℝ) →ₗ[ℝ] (ℝ → ℝ) →ₗ[ℝ] ℝ) (lc ls : (ℝ → ℝ) →ₗ[ℝ] ℝ),
      (∀ x y, B (tau x) (tau y) = K (x - y)) ∧ (∀ x, lc (tau x) = fc x) ∧
        ∀ x, ls (tau x) = fs x := by
  obtain ⟨B, hB⟩ := exists_bilinear_on_li tau_linearIndependent (fun x y => K (x - y))
  obtain ⟨lc, hlc⟩ := exists_linear_on_li tau_linearIndependent fc
  obtain ⟨ls, hls⟩ := exists_linear_on_li tau_linearIndependent fs
  exact ⟨B, lc, ls, hB, hlc, hls⟩

/-- **`WeilData` is non-empty.** The data are `Kψ` (pole, primes, Γ) and the cosh / sinh
values; no zero of ζ is an input. -/
theorem weilData_nonempty : Nonempty WeilData := by
  obtain ⟨B, lc, ls, hB, hlc, hls⟩ := exists_translate_data Kpsi
    (fun x => bigPsiI2 * Real.cosh (x / 2)) (fun x => bigPsiI2 * Real.sinh (x / 2))
  exact ⟨⟨B, lc, ls, hB, hlc, hls⟩⟩

/-- A chosen witness (CONSTRUCTED). -/
noncomputable def weilData : WeilData := Classical.choice weilData_nonempty

/-! ## 4. `HPsiEquiv` (PROVED, no premise) -/

theorem bigPsiI2_pos : 0 < bigPsiI2 := by
  unfold bigPsiI2
  have : 0 < Real.sinh (1 / 4) := Real.sinh_pos_iff.2 (by norm_num)
  positivity

/-- Every element of `F_ψ` is a finite combination `Σ cᵢ τ_{xᵢ}` indexed by `Fin n`. -/
theorem exists_fin_of_mem_Fpsi {v : ℝ → ℝ} (hv : v ∈ Fpsi) :
    ∃ (n : ℕ) (x c : Fin n → ℝ), v = ∑ i, c i • tau (x i) := by
  induction hv using Submodule.span_induction with
  | mem v hv =>
    obtain ⟨x, rfl⟩ := hv
    exact ⟨1, fun _ => x, fun _ => 1, by simp⟩
  | zero => exact ⟨0, Fin.elim0, Fin.elim0, by simp⟩
  | add v w _ _ hv hw =>
    obtain ⟨n, x, c, rfl⟩ := hv
    obtain ⟨m, y, d, rfl⟩ := hw
    exact ⟨n + m, Fin.append x y, Fin.append c d, by simp [Fin.sum_univ_add]⟩
  | smul a v _ hv =>
    obtain ⟨n, x, c, rfl⟩ := hv
    exact ⟨n, x, fun i => a * c i, by simp [Finset.smul_sum, smul_smul]⟩

theorem sum_mem_Fpsi {n : ℕ} (x c : Fin n → ℝ) : (∑ i, c i • tau (x i)) ∈ Fpsi :=
  Submodule.sum_mem _ fun i _ => Submodule.smul_mem _ _ (Submodule.subset_span ⟨x i, rfl⟩)

theorem B_sum (W : WeilData) {n : ℕ} (x c : Fin n → ℝ) :
    W.B (∑ i, c i • tau (x i)) (∑ i, c i • tau (x i)) = kform Kpsi x c := by
  simp only [map_sum, map_smul, LinearMap.sum_apply, LinearMap.smul_apply, W.hK,
    smul_eq_mul, kform, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
  ring

theorem lc_sum (W : WeilData) {n : ℕ} (x c : Fin n → ℝ) :
    W.lc (∑ i, c i • tau (x i)) = bigPsiI2 * ∑ i, c i * Real.cosh (x i / 2) := by
  simp only [map_sum, map_smul, W.hlc, smul_eq_mul, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  ring

/-- **The windowed form equals the Toeplitz form, for every `WeilData`** (PROVED).
`⇐`: a finite family fits in one window (`fit_translates`). `⇒`: an element of `F_ψ` is a
finite combination, `B` and `ℓ_c` are fixed on it by `hK`, `hlc`, and `Ψ(i/2) ≠ 0`. -/
theorem hPsiWindows_iff (W : WeilData) : HPsiWindows W ↔ HPsi := by
  constructor
  · intro H n x c hc
    obtain ⟨L, hL⟩ := fit_translates triPsi 1 (fun t ht => triPsi_eq_zero_of_one_le ht.le) n x
    have hvS : (∑ i, c i • tau (x i)) ∈ Spsi L :=
      Submodule.mem_inf.2 ⟨sum_mem_Fpsi x c,
        Submodule.sum_mem _ fun i _ => Submodule.smul_mem _ _ (hL i)⟩
    have hlc : W.lc (∑ i, c i • tau (x i)) = 0 := by rw [lc_sum, hc, mul_zero]
    have h := H L _ hvS hlc
    simp only at h
    rw [hlc, B_sum] at h
    linarith
  · intro H L v hv hlc
    simp only at hlc ⊢
    obtain ⟨n, x, c, rfl⟩ := exists_fin_of_mem_Fpsi (Submodule.mem_inf.1 hv).1
    rw [lc_sum] at hlc
    have hc : ∑ i, c i * Real.cosh (x i / 2) = 0 :=
      (mul_eq_zero.1 hlc).resolve_left bigPsiI2_pos.ne'
    have h := H n x c hc
    rw [lc_sum, hc, mul_zero, B_sum]
    linarith

/-- **`HPsiEquiv` holds** (PROVED). It is no longer a premise. -/
theorem hPsiEquiv : HPsiEquiv := hPsiWindows_iff

/-- **Non-vacuous form**: there is Weil data, and for it the windowed form is exactly `H_ψ`.
`H_ψ` itself remains OPEN. -/
theorem exists_weilData_iff : ∃ W : WeilData, (HPsiWindows W ↔ HPsi) :=
  ⟨weilData, hPsiWindows_iff weilData⟩

end GppKreinHPsi

#print axioms GppKreinHPsi.tau_linearIndependent
#print axioms GppKreinHPsi.exists_linear_on_li
#print axioms GppKreinHPsi.exists_bilinear_on_li
#print axioms GppKreinHPsi.exists_translate_data
#print axioms GppKreinHPsi.weilData_nonempty
#print axioms GppKreinHPsi.bigPsiI2_pos
#print axioms GppKreinHPsi.exists_fin_of_mem_Fpsi
#print axioms GppKreinHPsi.sum_mem_Fpsi
#print axioms GppKreinHPsi.B_sum
#print axioms GppKreinHPsi.lc_sum
#print axioms GppKreinHPsi.hPsiWindows_iff
#print axioms GppKreinHPsi.hPsiEquiv
#print axioms GppKreinHPsi.exists_weilData_iff
